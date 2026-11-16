"""Cross-file validation for a batch of individually-validated VCF files.

These checks can only run once every file has been parsed on its own:

1. reference version declarations must agree;
2. the same coordinate must carry the same REF in every file, otherwise the
   whole batch fails (no per-site fallback is allowed);
3. contigs must appear in the same order in every header, so a single
   coordinate sort can merge the streams deterministically.
"""
from __future__ import annotations

from dataclasses import dataclass

import pysam


@dataclass
class CrossCheckResult:
    valid: bool
    errors: list[str]
    reference_version: str | None


def cross_check(file_paths: list[tuple[str, str]]) -> CrossCheckResult:
    """``file_paths`` is a list of (display_name, stored_path)."""
    errors: list[str] = []
    reference_version: str | None = None

    refs_by_coord: dict[tuple[str, int], tuple[str, str]] = {}
    contig_orders: dict[str, list[str]] = {}

    opened: list[pysam.VariantFile] = []
    try:
        for name, path in file_paths:
            try:
                vcf = pysam.VariantFile(path)
            except (OSError, ValueError) as exc:
                errors.append(f"{name}：无法打开（{exc}）")
                continue
            opened.append(vcf)

            order = list(vcf.header.contigs.keys())
            contig_orders[name] = order

            declared = None
            for line_record in vcf.header.records:
                line = str(line_record)
                if line.startswith("##reference="):
                    declared = line[len("##reference="):].strip().strip('"')
                    break
            if declared:
                if reference_version is None:
                    reference_version = declared
                elif declared != reference_version:
                    errors.append(
                        f"参考版本不一致：{name} 声明 {declared!r}，"
                        f"此前文件声明 {reference_version!r}"
                    )

            try:
                for record in vcf:
                    key = (record.chrom, record.pos)
                    if key in refs_by_coord:
                        other_name, other_ref = refs_by_coord[key]
                        if other_ref != record.ref:
                            errors.append(
                                f"坐标 {record.chrom}:{record.pos} 的 REF 冲突："
                                f"{other_name} 为 {other_ref}，{name} 为 {record.ref}，"
                                "同坐标 REF 冲突导致整份合并失败"
                            )
                    else:
                        refs_by_coord[key] = (name, record.ref)
            except (OSError, ValueError) as exc:
                errors.append(f"{name}：遍历记录失败（{exc}）")
    finally:
        for vcf in opened:
            vcf.close()

    # Contig order consistency.
    baseline_name = next(iter(contig_orders), None)
    if baseline_name is not None:
        baseline = contig_orders[baseline_name]
        for name, order in contig_orders.items():
            if order != baseline:
                errors.append(
                    f"contig 顺序不一致：{baseline_name} 为 {baseline}，"
                    f"{name} 为 {order}，无法统一坐标排序"
                )
                break

    # De-duplicate while keeping the message order.
    seen: set[str] = set()
    unique_errors = [e for e in errors if not (e in seen or seen.add(e))]
    return CrossCheckResult(
        valid=not unique_errors,
        errors=unique_errors,
        reference_version=reference_version,
    )
