"""Workflow orchestration: uploads, validation, selections and merge workers.

Selection versioning is the core integrity mechanism:

* every batch has a monotonically increasing ``selection_version``;
* every generation records the version it was started with;
* any choice change bumps the version and immediately marks every existing
* generation -- completed or still running -- as ``superseded``;
* a worker checks the live version between sites and aborts, never publishing
  an output belonging to an outdated choice, even if the user changes back.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import threading
from pathlib import Path

from . import config, cross_check, db
from .merge import MergeAborted, SampleSource, merge_vcfs
from .validation import validate_vcf

logger = logging.getLogger(__name__)


class SelectionError(ValueError):
    """Raised when a generation is requested but prerequisites are unmet."""


# --------------------------------------------------------------------------- #
# Upload + validation
# --------------------------------------------------------------------------- #
def save_upload(batch_id: str, index: int, filename: str, fileobj) -> tuple[str, str]:
    lower = filename.lower()
    if lower.endswith(".vcf.gz"):
        suffix = ".vcf.gz"
    else:
        suffix = Path(filename).suffix.lower()
    stem = Path(filename).name[: -len(suffix)] if suffix else Path(filename).stem
    safe_stem = "".join(
        ch if ch.isalnum() or ch in "-_." else "_" for ch in stem
    )
    stored_name = f"{batch_id}_{index}_{safe_stem}{suffix}"
    stored_path = config.UPLOAD_DIR / stored_name
    config.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    with open(stored_path, "wb") as handle:
        shutil.copyfileobj(fileobj, handle)
    return filename, str(stored_path)


def process_uploaded_batch(batch_id: str, saved: list[tuple[str, str]]) -> None:
    """Validate all files, build duplicate groups and run cross-file checks."""
    with db.get_db() as conn:
        file_rows: list = []
        for original_name, stored_path in saved:
            file_id = db.create_file(conn, batch_id, original_name, stored_path)
            result = validate_vcf(stored_path)
            db.update_file_validation(
                conn,
                file_id,
                valid=result.valid,
                error=result.error,
                reference_version=result.reference_version,
                sample_count=len(result.samples),
                record_count=result.record_count,
                samples=result.samples,
                contigs=result.contigs,
            )
            file_rows.append(db.get_file(conn, file_id))

        valid_rows = [row for row in file_rows if row["valid"]]

        # Duplicate sample groups are derived from valid files only. Groups
        # with a single source are auto-selected: that original file is the
        # only possible choice. Multi-file groups start unresolved and block
        # generation until the operator picks one.
        sample_to_files: dict[str, list[str]] = {}
        for row in valid_rows:
            for sample in json.loads(row["samples_json"]):
                sample_to_files.setdefault(sample, []).append(row["id"])
        for sample, file_ids in sorted(sample_to_files.items()):
            group_id = db.create_dup_group(conn, batch_id, sample, file_ids)
            if len(file_ids) == 1:
                db.set_dup_group_choice(conn, group_id, file_ids[0])

        batch_errors: list[str] = []

        # Cross-file checks only over valid files.
        if len(valid_rows) >= 2:
            check = cross_check.cross_check(
                [(row["original_name"], row["stored_path"]) for row in valid_rows]
            )
            if check.reference_version:
                db.set_batch_reference(conn, batch_id, check.reference_version)
            batch_errors.extend(check.errors)
        else:
            invalid_names = [row["original_name"] for row in file_rows if not row["valid"]]
            batch_errors.append(
                "可参与合并的有效文件不足 2 份"
                + (f"，失效文件：{', '.join(invalid_names)}" if invalid_names else "")
            )

        if batch_errors:
            db.add_batch_errors(conn, batch_id, batch_errors)


# --------------------------------------------------------------------------- #
# Selection helpers
# --------------------------------------------------------------------------- #
def selection_summary(conn, batch_id: str) -> dict:
    groups = db.get_dup_groups(conn, batch_id)
    unresolved = [g["sample_name"] for g in groups if not g["chosen_file_id"]]
    return {
        "groups": groups,
        "unresolved": unresolved,
        "all_resolved": not unresolved,
    }


def assert_generation_allowed(conn, batch_id: str) -> None:
    batch = db.get_batch(conn, batch_id)
    if batch is None:
        raise SelectionError("批次不存在")
    if not batch["valid"]:
        errors = json.loads(batch["errors_json"])
        raise SelectionError(
            "批次存在阻断性错误，无法生成：" + "；".join(errors)
        )
    files = db.get_files(conn, batch_id)
    if sum(1 for f in files if f["valid"]) < config.MIN_FILES:
        raise SelectionError("有效 VCF 文件不足 2 份，无法生成")
    summary = selection_summary(conn, batch_id)
    if not summary["all_resolved"]:
        names = "、".join(summary["unresolved"])
        raise SelectionError(f"以下重复样本尚未选择唯一取用原件：{names}")
    active = db.get_active_generation(conn, batch_id)
    if active is not None:
        raise SelectionError("已有一个生成任务在进行中，请等待其结束")


def build_sources(conn, batch_id: str) -> list[SampleSource]:
    groups = db.get_dup_groups(conn, batch_id)
    sources: list[SampleSource] = []
    for group in groups:
        chosen = db.get_file(conn, group["chosen_file_id"])
        if chosen is None:
            raise SelectionError(
                f"样本 {group['sample_name']} 未选择取用原件"
            )
        sources.append(
            SampleSource(
                sample=group["sample_name"],
                file_path=chosen["stored_path"],
                file_name=chosen["original_name"],
            )
        )
    return sources


# --------------------------------------------------------------------------- #
# Background worker
# --------------------------------------------------------------------------- #
def _run_generation(gen_id: str, batch_id: str, version: int) -> None:
    with db.get_db() as conn:
        db.mark_generation(conn, gen_id, config.STATUS_RUNNING, started=True)
        try:
            sources = build_sources(conn, batch_id)
        except SelectionError as exc:
            db.mark_generation(
                conn, gen_id, config.STATUS_FAILED, error=str(exc), finished=True
            )
            return

    output_path = str(config.OUTPUT_DIR / f"cohort_{batch_id}_{gen_id}.vcf")

    def is_aborted() -> bool:
        with db.get_db() as check_conn:
            batch = db.get_batch(check_conn, batch_id)
            gen = db.get_generation(check_conn, gen_id)
            if batch is None or gen is None:
                return True
            if gen["status"] != config.STATUS_RUNNING:
                return True
            return batch["selection_version"] != version

    try:
        stats = merge_vcfs(sources, output_path, is_aborted=is_aborted)
    except MergeAborted:
        _safe_remove(output_path)
        with db.get_db() as conn:
            current = db.get_generation(conn, gen_id)
            if current is not None and current["status"] == config.STATUS_RUNNING:
                db.mark_generation(
                    conn,
                    gen_id,
                    config.STATUS_SUPERSEDED,
                    error="取用源选择已变更，该在途生成立即失效",
                    finished=True,
                )
        return
    except Exception as exc:  # noqa: BLE001 - failure reason must be surfaced
        _safe_remove(output_path)
        with db.get_db() as conn:
            db.mark_generation(
                conn,
                gen_id,
                config.STATUS_FAILED,
                error=f"合并失败：{exc}",
                finished=True,
            )
        return

    # Re-check before publishing: never publish a stale result.
    if is_aborted():
        _safe_remove(output_path)
        with db.get_db() as conn:
            current = db.get_generation(conn, gen_id)
            if current is not None and current["status"] == config.STATUS_RUNNING:
                db.mark_generation(
                    conn,
                    gen_id,
                    config.STATUS_SUPERSEDED,
                    error="取用源选择已变更，结果不予发布",
                    finished=True,
                )
        return

    logger.info(
        "generation %s merged %s sites for %s samples -> %s",
        gen_id, stats["sites"], stats["samples"], output_path,
    )
    with db.get_db() as conn:
        db.mark_generation(
            conn,
            gen_id,
            config.STATUS_DONE,
            output_path=output_path,
            error="",
            finished=True,
        )


def start_generation(batch_id: str) -> str:
    """Validate prerequisites synchronously, then spawn the worker thread."""
    with db.get_db() as conn:
        assert_generation_allowed(conn, batch_id)
        version = db.get_batch(conn, batch_id)["selection_version"]
        gen_id = db.create_generation(conn, batch_id, version)

    thread = threading.Thread(
        target=_run_generation,
        args=(gen_id, batch_id, version),
        name=f"merge-{gen_id[:8]}",
        daemon=True,
    )
    thread.start()
    return gen_id


def change_selection(batch_id: str, group_id: str, file_id: str | None) -> None:
    """Change a sample's source; invalidates all generations immediately."""
    with db.get_db() as conn:
        group = db.get_dup_group(conn, group_id)
        if group is None or group["batch_id"] != batch_id:
            raise SelectionError("重复样本分组不存在")
        updated = db.set_dup_group_choice(conn, group_id, file_id)
        if updated is None:
            raise SelectionError("选择无效")
        db.bump_selection_version(conn, batch_id)

    # Remove output files of generations that were just superseded so no
    # half-product or stale result can ever be downloaded.
    _cleanup_invalid_outputs(batch_id)


def _cleanup_invalid_outputs(batch_id: str) -> None:
    with db.get_db() as conn:
        rows = conn.execute(
            "SELECT output_path FROM generations WHERE batch_id=? "
            "AND status=? AND output_path IS NOT NULL",
            (batch_id, config.STATUS_SUPERSEDED),
        ).fetchall()
    for row in rows:
        _safe_remove(row["output_path"])


def _safe_remove(path: str | None) -> None:
    if not path:
        return
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
