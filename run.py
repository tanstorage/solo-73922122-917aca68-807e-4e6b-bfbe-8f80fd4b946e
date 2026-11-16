"""Local development entry point: ``python run.py``."""
import uvicorn

from app.config import ensure_dirs
from app.db import init_db

if __name__ == "__main__":
    ensure_dirs()
    init_db()
    uvicorn.run("app.main:app", host="0.0.0.0", port=8000, reload=False)
