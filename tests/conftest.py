"""Pytest fixtures: isolated data directory and FastAPI test client.

The application reads ``VCF_DATA_DIR`` lazily through :mod:`app.config`, so
the environment variable is set before the app package is imported. Modules
are intentionally never reloaded: reloading would give request handlers and
background worker threads different module objects (and therefore different
thread-local connection registries), which breaks end-to-end tests.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


@pytest.fixture()
def app_env(tmp_path, monkeypatch):
    """Point the application at a throwaway data directory."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    monkeypatch.setenv("VCF_DATA_DIR", str(data_dir))

    # Import after the env var is set; config paths are evaluated at import.
    from app import config, db, main, merge, service, validation

    config.ensure_dirs()
    db.init_db()
    yield {
        "config": config,
        "db": db,
        "service": service,
        "main": main,
        "merge": merge,
        "validation": validation,
        "data_dir": data_dir,
    }


@pytest.fixture()
def client(app_env):
    from fastapi.testclient import TestClient

    with TestClient(app_env["main"].app) as test_client:
        yield test_client


DEMO_DIR = ROOT / "demo"
