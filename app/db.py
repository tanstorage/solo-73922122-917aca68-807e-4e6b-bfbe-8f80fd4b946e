"""SQLite persistence layer.

Schema overview
---------------
batches            : one upload transaction (2-3 VCF files)
vcf_files          : uploaded files and their per-file validation outcome
dup_groups         : groups of files sharing a sample; exactly one must be
                     chosen as the unique source of that sample
generations        : merge runs, tagged with the ``selection_version`` they
                     were started with; any later selection change invalidates
                     them immediately

The schema deliberately keeps no per-locus data: VCF contents stay on disk
and the database only tracks workflow state.
"""
from __future__ import annotations

import json
import sqlite3
import threading
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from . import config

_local = threading.local()
_init_lock = threading.Lock()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def new_id() -> str:
    return uuid.uuid4().hex


def _connect(db_path: Path | None = None) -> sqlite3.Connection:
    conn = sqlite3.connect(
        str(db_path or config.DB_PATH),
        timeout=30,
        check_same_thread=False,
        isolation_level=None,  # autocommit; explicit BEGIN for transactions
    )
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA busy_timeout=30000")
    conn.execute("PRAGMA synchronous=NORMAL")
    return conn


@contextmanager
def get_db() -> Iterator[sqlite3.Connection]:
    """Thread-local connection used by request handlers and workers."""
    conn = getattr(_local, "conn", None)
    if conn is None:
        conn = _connect()
        _local.conn = conn
    try:
        yield conn
    except sqlite3.Error:
        conn.rollback()
        raise


SCHEMA = """
CREATE TABLE IF NOT EXISTS batches (
    id                 TEXT PRIMARY KEY,
    created_at         TEXT NOT NULL,
    reference_version  TEXT,
    valid              INTEGER NOT NULL DEFAULT 1,
    -- Batch-level blocking reasons (reference mismatch, REF conflicts ...)
    errors_json        TEXT NOT NULL DEFAULT '[]',
    selection_version  INTEGER NOT NULL DEFAULT 1
);

CREATE TABLE IF NOT EXISTS vcf_files (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    original_name      TEXT NOT NULL,
    stored_path        TEXT NOT NULL,
    uploaded_at        TEXT NOT NULL,
    valid              INTEGER NOT NULL DEFAULT 1,
    error              TEXT,
    reference_version  TEXT,
    sample_count       INTEGER NOT NULL DEFAULT 0,
    record_count       INTEGER NOT NULL DEFAULT 0,
    samples_json       TEXT NOT NULL DEFAULT '[]',
    contigs_json       TEXT NOT NULL DEFAULT '[]'
);

CREATE TABLE IF NOT EXISTS dup_groups (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    sample_name        TEXT NOT NULL,
    -- Files that carry this sample (vcf_files id list, JSON array)
    file_ids_json      TEXT NOT NULL,
    -- Chosen source file id; NULL means no selection yet (blocks generation)
    chosen_file_id     TEXT
);

CREATE TABLE IF NOT EXISTS generations (
    id                 TEXT PRIMARY KEY,
    batch_id           TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    status             TEXT NOT NULL,
    selection_version  INTEGER NOT NULL,
    output_path        TEXT,
    error              TEXT,
    created_at         TEXT NOT NULL,
    started_at         TEXT,
    finished_at        TEXT
);

CREATE INDEX IF NOT EXISTS idx_files_batch ON vcf_files(batch_id);
CREATE INDEX IF NOT EXISTS idx_dup_batch ON dup_groups(batch_id);
CREATE INDEX IF NOT EXISTS idx_gen_batch ON generations(batch_id);
"""


def init_db() -> None:
    with _init_lock:
        config.ensure_dirs()
        conn = _connect()
        try:
            conn.executescript(SCHEMA)
            conn.execute(
                "UPDATE generations SET status=?, finished_at=?, "
                "error=COALESCE(error, ?) WHERE status IN (?, ?)",
                (
                    config.STATUS_FAILED,
                    _now(),
                    "服务在合并完成前重启，该在途生成判定失败",
                    config.STATUS_PENDING,
                    config.STATUS_RUNNING,
                ),
            )
        finally:
            conn.close()


# --------------------------------------------------------------------------- #
# Reads
# --------------------------------------------------------------------------- #
def get_batch(conn: sqlite3.Connection, batch_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM batches WHERE id=?", (batch_id,)
    ).fetchone()


def get_files(conn: sqlite3.Connection, batch_id: str) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM vcf_files WHERE batch_id=? ORDER BY uploaded_at, id",
        (batch_id,),
    ).fetchall()


def get_file(conn: sqlite3.Connection, file_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM vcf_files WHERE id=?", (file_id,)
    ).fetchone()


def get_dup_groups(conn: sqlite3.Connection, batch_id: str) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM dup_groups WHERE batch_id=? ORDER BY sample_name",
        (batch_id,),
    ).fetchall()


def get_dup_group(conn: sqlite3.Connection, group_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM dup_groups WHERE id=?", (group_id,)
    ).fetchone()


def get_active_generation(
    conn: sqlite3.Connection, batch_id: str
) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM generations WHERE batch_id=? AND status IN (?, ?) "
        "ORDER BY created_at DESC LIMIT 1",
        (batch_id, config.STATUS_PENDING, config.STATUS_RUNNING),
    ).fetchone()


def get_latest_generation(
    conn: sqlite3.Connection, batch_id: str
) -> sqlite3.Row | None:
    # rowid reflects insertion order; created_at has only second-level
    # resolution and several generations can share one timestamp.
    return conn.execute(
        "SELECT * FROM generations WHERE batch_id=? ORDER BY rowid DESC "
        "LIMIT 1",
        (batch_id,),
    ).fetchone()


def get_generation(conn: sqlite3.Connection, gen_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM generations WHERE id=?", (gen_id,)
    ).fetchone()


def get_running_generations(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM generations WHERE status IN (?, ?)",
        (config.STATUS_PENDING, config.STATUS_RUNNING),
    ).fetchall()


# --------------------------------------------------------------------------- #
# Writes
# --------------------------------------------------------------------------- #
def create_batch(conn: sqlite3.Connection) -> str:
    batch_id = new_id()
    conn.execute(
        "INSERT INTO batches(id, created_at) VALUES(?, ?)",
        (batch_id, _now()),
    )
    return batch_id


def create_file(
    conn: sqlite3.Connection,
    batch_id: str,
    original_name: str,
    stored_path: str,
) -> str:
    file_id = new_id()
    conn.execute(
        "INSERT INTO vcf_files(id, batch_id, original_name, stored_path, "
        "uploaded_at) VALUES(?, ?, ?, ?, ?)",
        (file_id, batch_id, original_name, stored_path, _now()),
    )
    return file_id


def update_file_validation(
    conn: sqlite3.Connection,
    file_id: str,
    *,
    valid: bool,
    error: str | None,
    reference_version: str | None,
    sample_count: int,
    record_count: int,
    samples: list[str],
    contigs: list[str],
) -> None:
    conn.execute(
        "UPDATE vcf_files SET valid=?, error=?, reference_version=?, "
        "sample_count=?, record_count=?, samples_json=?, contigs_json=? "
        "WHERE id=?",
        (
            1 if valid else 0,
            error,
            reference_version,
            sample_count,
            record_count,
            json.dumps(samples),
            json.dumps(contigs),
            file_id,
        ),
    )


def set_batch_reference(conn: sqlite3.Connection, batch_id: str, ref: str | None) -> None:
    conn.execute(
        "UPDATE batches SET reference_version=? WHERE id=?", (ref, batch_id)
    )


def add_batch_errors(conn: sqlite3.Connection, batch_id: str, errors: list[str]) -> None:
    row = get_batch(conn, batch_id)
    existing: list[str] = json.loads(row["errors_json"]) if row else []
    existing.extend(errors)
    conn.execute(
        "UPDATE batches SET errors_json=?, valid=? WHERE id=?",
        (json.dumps(existing), 0 if existing else 1, batch_id),
    )


def create_dup_group(
    conn: sqlite3.Connection,
    batch_id: str,
    sample_name: str,
    file_ids: list[str],
) -> str:
    group_id = new_id()
    conn.execute(
        "INSERT INTO dup_groups(id, batch_id, sample_name, file_ids_json, "
        "chosen_file_id) VALUES(?, ?, ?, ?, NULL)",
        (group_id, batch_id, sample_name, json.dumps(file_ids)),
    )
    return group_id


def bump_selection_version(conn: sqlite3.Connection, batch_id: str) -> int:
    """Invalidate every result/in-flight run tied to older selections."""
    conn.execute("BEGIN IMMEDIATE")
    row = get_batch(conn, batch_id)
    new_version = row["selection_version"] + 1
    conn.execute(
        "UPDATE batches SET selection_version=? WHERE id=?",
        (new_version, batch_id),
    )
    conn.execute(
        "UPDATE generations SET status=?, finished_at=?, error=? "
        "WHERE batch_id=? AND status IN (?, ?)",
        (
            config.STATUS_SUPERSEDED,
            _now(),
            "取用源选择已变更，旧结果与在途生成立即失效",
            batch_id,
            config.STATUS_PENDING,
            config.STATUS_RUNNING,
        ),
    )
    # Completed outputs of previous selections become invalid as well.
    conn.execute(
        "UPDATE generations SET status=?, finished_at=?, error=? "
        "WHERE batch_id=? AND status=?",
        (
            config.STATUS_SUPERSEDED,
            _now(),
            "取用源选择已变更，旧结果与在途生成立即失效",
            batch_id,
            config.STATUS_DONE,
        ),
    )
    conn.execute("COMMIT")
    return new_version


def set_dup_group_choice(
    conn: sqlite3.Connection, group_id: str, file_id: str | None
) -> sqlite3.Row | None:
    """Persist a chosen source; returns the updated row or None on bad input."""
    group = get_dup_group(conn, group_id)
    if group is None:
        return None
    if file_id is not None:
        allowed = json.loads(group["file_ids_json"])
        if file_id not in allowed:
            raise ValueError("chosen file is not an original source of this sample")
    conn.execute(
        "UPDATE dup_groups SET chosen_file_id=? WHERE id=?",
        (file_id, group_id),
    )
    return get_dup_group(conn, group_id)


def create_generation(
    conn: sqlite3.Connection, batch_id: str, selection_version: int
) -> str:
    gen_id = new_id()
    conn.execute(
        "INSERT INTO generations(id, batch_id, status, selection_version, "
        "created_at) VALUES(?, ?, ?, ?, ?)",
        (gen_id, batch_id, config.STATUS_PENDING, selection_version, _now()),
    )
    return gen_id


def mark_generation(
    conn: sqlite3.Connection,
    gen_id: str,
    status: str,
    *,
    output_path: str | None = None,
    error: str | None = None,
    started: bool = False,
    finished: bool = False,
) -> None:
    sets = ["status=?"]
    params: list[Any] = [status]
    if output_path is not None:
        sets.append("output_path=?")
        params.append(output_path)
    if error is not None:
        sets.append("error=?")
        params.append(error)
    if started:
        sets.append("started_at=?")
        params.append(_now())
    if finished:
        sets.append("finished_at=?")
        params.append(_now())
    params.append(gen_id)
    conn.execute(f"UPDATE generations SET {', '.join(sets)} WHERE id=?", params)
