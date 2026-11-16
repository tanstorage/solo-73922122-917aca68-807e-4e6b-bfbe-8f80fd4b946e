"""FastAPI application: VCF cohort merge web UI and JSON API."""
from __future__ import annotations

import json
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles

from . import config, db, service
from .service import SelectionError


@asynccontextmanager
async def lifespan(_app: FastAPI):
    config.ensure_dirs()
    db.init_db()
    yield


app = FastAPI(
    title="VCF 队列合并服务",
    description="分批检测二倍体 VCF 的队列级合并，供关联分析使用",
    version="1.0.0",
    lifespan=lifespan,
)

STATIC_DIR = Path(__file__).resolve().parent / "static"


# --------------------------------------------------------------------------- #
# Serialization helpers
# --------------------------------------------------------------------------- #
def _file_payload(row) -> dict:
    return {
        "id": row["id"],
        "original_name": row["original_name"],
        "valid": bool(row["valid"]),
        "error": row["error"],
        "reference_version": row["reference_version"],
        "sample_count": row["sample_count"],
        "record_count": row["record_count"],
        "samples": json.loads(row["samples_json"]),
        "contigs": json.loads(row["contigs_json"]),
        "uploaded_at": row["uploaded_at"],
    }


def _group_payload(conn, row) -> dict:
    file_ids = json.loads(row["file_ids_json"])
    choices = []
    for file_id in file_ids:
        file_row = db.get_file(conn, file_id)
        choices.append(
            {
                "file_id": file_id,
                "file_name": file_row["original_name"] if file_row else file_id,
            }
        )
    chosen = row["chosen_file_id"]
    return {
        "id": row["id"],
        "sample_name": row["sample_name"],
        "choices": choices,
        "chosen_file_id": chosen,
        "is_duplicate": len(file_ids) > 1,
    }


def _generation_payload(row) -> dict:
    return {
        "id": row["id"],
        "status": row["status"],
        "selection_version": row["selection_version"],
        "error": row["error"],
        "output_path": row["output_path"],
        "created_at": row["created_at"],
        "started_at": row["started_at"],
        "finished_at": row["finished_at"],
        "downloadable": bool(
            row["status"] == config.STATUS_DONE
            and row["output_path"]
            and Path(row["output_path"]).is_file()
        ),
    }


def _batch_payload(conn, batch_id: str) -> dict:
    batch = db.get_batch(conn, batch_id)
    if batch is None:
        raise HTTPException(status_code=404, detail="批次不存在")
    files = [_file_payload(row) for row in db.get_files(conn, batch_id)]
    groups = [_group_payload(conn, row) for row in db.get_dup_groups(conn, batch_id)]
    latest = db.get_latest_generation(conn, batch_id)
    unresolved = [g["sample_name"] for g in groups if not g["chosen_file_id"]]
    return {
        "id": batch["id"],
        "created_at": batch["created_at"],
        "reference_version": batch["reference_version"],
        "batch_valid": bool(batch["valid"]),
        "batch_errors": json.loads(batch["errors_json"]),
        "selection_version": batch["selection_version"],
        "files": files,
        "dup_groups": groups,
        "unresolved_samples": unresolved,
        "ready_to_generate": bool(batch["valid"]) and not unresolved
        and sum(1 for f in files if f["valid"]) >= config.MIN_FILES,
        "latest_generation": _generation_payload(latest) if latest else None,
    }


# --------------------------------------------------------------------------- #
# API endpoints
# --------------------------------------------------------------------------- #
@app.post("/api/batches")
async def create_batch(files: list[UploadFile] = File(...)):
    if not config.MIN_FILES <= len(files) <= config.MAX_FILES:
        raise HTTPException(
            status_code=400,
            detail=f"每次必须上传 {config.MIN_FILES}~{config.MAX_FILES} 份 VCF 文件",
        )

    saved: list[tuple[str, str]] = []
    with db.get_db() as conn:
        batch_id = db.create_batch(conn)

    try:
        for index, upload in enumerate(files):
            filename = upload.filename or ""
            if not any(
                filename.lower().endswith(ext)
                for ext in config.ALLOWED_EXTENSIONS
            ):
                raise HTTPException(
                    status_code=400,
                    detail=(
                        f"文件 {filename} 扩展名不支持，"
                        "仅允许 .vcf / .vcf.gz"
                    ),
                )
            # Guard against oversized uploads.
            config.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
            original_name, stored_path = service.save_upload(
                batch_id, index, upload.filename, upload.file
            )
            size = Path(stored_path).stat().st_size
            if size > config.MAX_UPLOAD_BYTES:
                Path(stored_path).unlink(missing_ok=True)
                raise HTTPException(
                    status_code=400,
                    detail=f"文件 {upload.filename} 超过大小限制",
                )
            saved.append((original_name, stored_path))

        service.process_uploaded_batch(batch_id, saved)
    except HTTPException:
        # Roll back partial uploads on rejection.
        for _, stored_path in saved:
            Path(stored_path).unlink(missing_ok=True)
        with db.get_db() as conn:
            conn.execute("DELETE FROM batches WHERE id=?", (batch_id,))
        raise
    except Exception as exc:  # noqa: BLE001
        # Keep the batch but record a blocking reason so the UI can show it
        # instead of leaving an orphan batch with no explanation.
        with db.get_db() as conn:
            db.add_batch_errors(
                conn,
                batch_id,
                [f"上传处理过程中发生未预期错误：{exc}"],
            )
        raise HTTPException(status_code=500, detail="服务器处理上传时发生内部错误")

    with db.get_db() as conn:
        return _batch_payload(conn, batch_id)


@app.get("/api/batches/{batch_id}")
def get_batch(batch_id: str):
    with db.get_db() as conn:
        return _batch_payload(conn, batch_id)


@app.put("/api/batches/{batch_id}/selection/{group_id}")
def update_selection(
    batch_id: str, group_id: str, file_id: str | None = Form(None)
):
    if file_id == "":
        file_id = None
    try:
        service.change_selection(batch_id, group_id, file_id)
    except SelectionError as exc:
        raise HTTPException(status_code=409, detail=str(exc))
    with db.get_db() as conn:
        return _batch_payload(conn, batch_id)


@app.post("/api/batches/{batch_id}/generate")
def generate(batch_id: str):
    try:
        service.start_generation(batch_id)
    except SelectionError as exc:
        raise HTTPException(status_code=409, detail=str(exc))
    with db.get_db() as conn:
        return _batch_payload(conn, batch_id)


@app.get("/api/batches/{batch_id}/download")
def download(batch_id: str):
    with db.get_db() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            raise HTTPException(status_code=404, detail="批次不存在")
        latest = db.get_latest_generation(conn, batch_id)
        if latest is None:
            raise HTTPException(status_code=409, detail="尚无生成结果")
        if latest["status"] != config.STATUS_DONE:
            raise HTTPException(
                status_code=409,
                detail="仅当前选择对应的完整 VCF 可下载，该结果已失效或未完成",
            )
        output_path = latest["output_path"]
        if not output_path or not Path(output_path).is_file():
            raise HTTPException(
                status_code=409, detail="结果文件不存在（可能因服务重启未完成）"
            )
        if latest["selection_version"] != batch["selection_version"]:
            raise HTTPException(
                status_code=409, detail="取用源已变更，该结果不再属于当前选择"
            )
    download_name = f"cohort_{batch_id}.vcf"
    return FileResponse(
        output_path,
        media_type="text/vcard",
        filename=download_name,
    )


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/", response_class=HTMLResponse)
def index():
    return FileResponse(str(STATIC_DIR / "index.html"))


app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
