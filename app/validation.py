"""Per-file VCF validation.

Rules enforced here (one failure marks the whole file invalid):

* the file must be parseable by pysam and declare at least one sample/site;
* every contig used by a record must be declared in the header and records
  must be coordinate-sorted following the header contig order;
* each coordinate may appear exactly once (no per-site splicing later);
* REF and every ALT must be non-empty ACGT sequences (ALT ``.`` means the
  record has no alternate allele and is accepted);
* FORMAT fields are restricted to GT/AD/PL and GT is mandatory;
* every sample genotype must be diploid (two allele slots, missing allowed).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pysam

from . import config


@dataclass
class FileValidation:
    valid: bool
    error: str | None = None
    reference_version: str | None = None
    samples: list[str] = field(default_factory=list)
    contigs: list[str] = field(default_factory=list)
    record_count: int = 0


def parse_reference_version(header: pysam.VariantHeader) -> str | None:
    """Extract the reference declaration from ``##reference=`` header lines."""
    for record in header.records:
        line = str(record)
        if line.startswith("##reference="):
            value = line[len("##reference="):].strip()
            # Values may be quoted or a file:// URL.
            value = value.strip('"')
            if value:
                return value
    return None


def _allele_ok(allele: str | None) -> bool:
    return allele is not None and config.ALLELE_RE.fullmatch(allele) is not None


def validate_vcf(path: str) -> FileValidation:
    try:
        vcf = pysam.VariantFile(path)
    except (OSError, ValueError, NotImplementedError) as exc:
        hint = ""
        if path.lower().endswith(".gz") and "seek" in str(exc).lower():
            hint = "（.vcf.gz 需为 bgzip 压缩的 BGZF 格式，普通 gzip 不被支持）"
        return FileValidation(valid=False, error=f"无法解析 VCF 文件{hint}：{exc}")

    try:
        header = vcf.header
        samples = list(header.samples)
        contigs = list(header.contigs.keys())
        reference_version = parse_reference_version(header)

        if not samples:
            return FileValidation(
                valid=False,
                error="文件不包含任何样本列（#CHROM 行缺少样本）",
                reference_version=reference_version,
                samples=samples,
                contigs=contigs,
            )

        contig_index = {name: idx for idx, name in enumerate(contigs)}
        seen_coords: set[tuple[str, int]] = set()
        prev_key: tuple[int, int] | None = None
        record_count = 0

        for record in vcf:
            coord = f"{record.chrom}:{record.pos}"

            if record.chrom not in contig_index:
                return FileValidation(
                    valid=False,
                    error=(
                        f"记录 {coord} 使用的 contig 未在文件头 ##contig 中声明，"
                        "无法保证坐标排序"
                    ),
                    reference_version=reference_version,
                    samples=samples,
                    contigs=contigs,
                    record_count=record_count,
                )

            key = (contig_index[record.chrom], record.pos)
            if prev_key is not None and key < prev_key:
                return FileValidation(
                    valid=False,
                    error=(
                        f"记录 {coord} 未按文件头 contig 顺序与坐标升序排列，"
                        "要求上传坐标排序的 VCF"
                    ),
                    reference_version=reference_version,
                    samples=samples,
                    contigs=contigs,
                    record_count=record_count,
                )
            prev_key = key

            if (record.chrom, record.pos) in seen_coords:
                return FileValidation(
                    valid=False,
                    error=f"坐标 {coord} 存在多条记录，要求每份文件同坐标仅一条记录",
                    reference_version=reference_version,
                    samples=samples,
                    contigs=contigs,
                    record_count=record_count,
                )
            seen_coords.add((record.chrom, record.pos))

            if not _allele_ok(record.ref):
                return FileValidation(
                    valid=False,
                    error=(
                        f"记录 {coord} 的 REF={record.ref!r} 不是非空 ACGT 序列，"
                        "等位基因仅限 ACGT"
                    ),
                    reference_version=reference_version,
                    samples=samples,
                    contigs=contigs,
                    record_count=record_count,
                )

            for alt in record.alts or ():
                if alt != "." and not _allele_ok(alt):
                    return FileValidation(
                        valid=False,
                        error=(
                            f"记录 {coord} 的 ALT={alt!r} 不是非空 ACGT 序列，"
                            "等位基因仅限 ACGT"
                        ),
                        reference_version=reference_version,
                        samples=samples,
                        contigs=contigs,
                        record_count=record_count,
                    )

            format_fields = list(record.format)
            for field_name in format_fields:
                if field_name not in config.ALLOWED_FORMAT_FIELDS:
                    return FileValidation(
                        valid=False,
                        error=(
                            f"记录 {coord} 的 FORMAT 含不支持字段 {field_name}，"
                            "仅允许 GT、AD、PL"
                        ),
                        reference_version=reference_version,
                        samples=samples,
                        contigs=contigs,
                        record_count=record_count,
                    )
            if "GT" not in format_fields:
                return FileValidation(
                    valid=False,
                    error=f"记录 {coord} 缺少必需的 GT 字段",
                    reference_version=reference_version,
                    samples=samples,
                    contigs=contigs,
                    record_count=record_count,
                )

            for sample_name in samples:
                gt = record.samples[sample_name].get("GT")
                if gt is not None and len(gt) != config.DIPLOID_PLOIDY:
                    ploidy = 0 if gt is None else len(gt)
                    return FileValidation(
                        valid=False,
                        error=(
                            f"样本 {sample_name} 在 {coord} 的倍性为 {ploidy}，"
                            "要求二倍体（允许 ./. 缺失）"
                        ),
                        reference_version=reference_version,
                        samples=samples,
                        contigs=contigs,
                        record_count=record_count,
                    )

            record_count += 1

        if record_count == 0:
            return FileValidation(
                valid=False,
                error="文件不含任何位点记录",
                reference_version=reference_version,
                samples=samples,
                contigs=contigs,
            )

        return FileValidation(
            valid=True,
            reference_version=reference_version,
            samples=samples,
            contigs=contigs,
            record_count=record_count,
        )
    except (OSError, ValueError) as exc:
        return FileValidation(
            valid=False,
            error=f"读取 VCF 记录时失败：{exc}",
            reference_version=reference_version if "reference_version" in locals() else None,
            samples=list(header.samples) if "header" in locals() else [],
            contigs=list(header.contigs.keys()) if "header" in locals() else [],
            record_count=record_count if "record_count" in locals() else 0,
        )
    finally:
        vcf.close()


def references_conflict(refs: list[tuple[str, str | None]]) -> str | None:
    """Return an error string if declared reference versions disagree.

    ``refs`` is a list of (file_name, reference_version). Files without a
    reference declaration are ignored; two explicit different values conflict.
    """
    declared = [(name, ref) for name, ref in refs if ref]
    distinct = {ref for _, ref in declared}
    if len(distinct) > 1:
        detail = "，".join(f"{name}={ref!r}" for name, ref in declared)
        return f"参考版本不一致（{detail}），要求同参考版本的 VCF"
    return None
