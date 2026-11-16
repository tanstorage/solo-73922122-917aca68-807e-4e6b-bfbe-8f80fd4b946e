"""Minimal BGZF (Blocked GZIF, RFC 1952) writer used for demo/test fixtures.

pysam/htslib requires ``.vcf.gz`` files to be BGZF rather than plain gzip,
because BGZF supports virtual offsets and seeking. This module produces
valid BGZF blocks using only the Python standard library (zlib), so no
external ``bgzip`` binary is needed.

The production application does not need this: users are expected to upload
files produced by ``bgzip``. It exists so the repository's demo and tests can
generate a conforming ``.vcf.gz`` without any non-Python toolchain.
"""
from __future__ import annotations

import struct
import zlib

# BGZF blocks are capped at 64 KiB including headers/trailer.
MAX_BLOCK_DATA = 65280


def _bgzf_block(data: bytes) -> bytes:
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    deflated = compressor.compress(data) + compressor.flush()
    # BSIZE = total block size minus 1.
    bsize = len(deflated) + 25
    header = b"".join(
        [
            b"\x1f\x8b",  # gzip magic
            b"\x08",  # compression method: deflate
            b"\x04",  # flags: FEXTRA
            b"\x00\x00\x00\x00",  # mtime
            b"\x00",  # extra flags
            b"\xff",  # OS: unknown
            struct.pack("<H", 6),  # extra field length
            b"BC",  # extra field IDs
            struct.pack("<H", 2),  # extra field size
            struct.pack("<H", bsize),  # BSIZE
        ]
    )
    crc = zlib.crc32(data) & 0xFFFFFFFF
    trailer = struct.pack("<II", crc, len(data))
    return header + deflated + trailer


def bgzip(src_path: str, dst_path: str) -> None:
    """Compress ``src_path`` into a BGZF file at ``dst_path``."""
    with open(src_path, "rb") as handle:
        raw = handle.read()
    with open(dst_path, "wb") as out:
        for offset in range(0, len(raw), MAX_BLOCK_DATA):
            out.write(_bgzf_block(raw[offset : offset + MAX_BLOCK_DATA]))
        # Mandatory empty BGZF end-of-stream marker.
        out.write(_bgzf_block(b""))


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        raise SystemExit("usage: python bgzf_util.py src.vcf dst.vcf.gz")
    bgzip(sys.argv[1], sys.argv[2])
