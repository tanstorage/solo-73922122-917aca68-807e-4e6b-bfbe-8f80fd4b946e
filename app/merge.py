"""Cohort VCF merge engine.

Given a set of source VCF files and an explicit per-sample source selection,
produce one cohort VCF where:

* every distinct coordinate appears exactly once with the union of ALT
  alleles (ordered by first appearance across files);
* GT allele indices are remapped onto the union;
* AD (one value per allele) and PL (one value per unordered diploid genotype)
  are remapped onto the union, values for alleles introduced by other files
  stay missing (``.``) instead of being fabricated;
* samples that do not have a site in their source file are emitted as fully
  missing (``./.``) at that site -- never silently converted to 0/0;
* the merge never splices sites: each sample's column comes from one chosen
  original file for the entire cohort file.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import pysam


class MergeAborted(Exception):
    """Raised when the worker notices its generation was invalidated."""


@dataclass(frozen=True)
class SampleSource:
    sample: str
    file_path: str
    file_name: str


def _pl_index(i: int, j: int, n_alleles: int) -> int:
    """VCF PL ordering for unordered diploid genotypes.

    Layout is (0,0),(0,1),(1,1),(0,2),(1,2),(2,2),(0,3)... so for i <= j the
    index is j*(j+1)/2 + i.
    """
    if i > j:
        i, j = j, i
    return j * (j + 1) // 2 + i


def _triangular_count(n_alleles: int) -> int:
    return n_alleles * (n_alleles + 1) // 2


@dataclass
class _Cursor:
    name: str
    path: str
    vcf: pysam.VariantFile
    record: pysam.VariantRecord | None
    contig_index: dict[str, int]

    def advance(self) -> None:
        try:
            self.record = next(self.vcf)
        except StopIteration:
            self.record = None

    def sort_key(self) -> tuple[int, int]:
        record = self.record
        assert record is not None
        return self.contig_index[record.chrom], record.pos


def _remap_genotype(record: pysam.VariantRecord, sample: str, allele_map: list[int | None]) -> tuple:
    """Remap one sample's GT onto the merged allele numbering."""
    sample_data = record.samples[sample]
    gt = sample_data.get("GT")
    if gt is None or None in gt:
        return (None, None)
    remapped = tuple(allele_map[idx] for idx in gt)
    if None in remapped:
        # Defensive: every source allele is part of the union, so this
        # should never happen for a valid input.
        return (None, None)
    return remapped  # type: ignore[return-value]


def _remap_ad(record: pysam.VariantRecord, sample: str, allele_map: list[int | None], n_out: int) -> tuple:
    values = record.samples[sample].get("AD")
    out: list[int | None] = [None] * n_out
    if values is None:
        return tuple(out)
    for src_idx, value in enumerate(values):
        if src_idx >= len(allele_map):
            # Malformed AD with more values than declared alleles: ignore
            # the surplus rather than index outside the mapping.
            break
        out_idx = allele_map[src_idx]
        if out_idx is not None and value is not None:
            out[out_idx] = value
    return tuple(out)


def _remap_pl(
    record: pysam.VariantRecord,
    sample: str,
    allele_map: list[int | None],
    n_out: int,
) -> tuple:
    values = record.samples[sample].get("PL")
    out: list[int | None] = [None] * _triangular_count(n_out)
    if values is None:
        return tuple(out)
    n_src = len(allele_map)
    for src_i in range(n_src):
        out_i = allele_map[src_i]
        if out_i is None:
            continue
        for src_j in range(src_i, n_src):
            out_j = allele_map[src_j]
            if out_j is None:
                continue
            src_index = _pl_index(src_i, src_j, n_src)
            if src_index >= len(values):
                break
            value = values[src_index]
            if value is None:
                continue
            out[_pl_index(out_i, out_j, n_out)] = value
    return tuple(out)


def merge_vcfs(
    sources: list[SampleSource],
    output_path: str,
    *,
    is_aborted=None,
) -> dict:
    """Merge selected samples from source files into a single cohort VCF.

    ``is_aborted`` is an optional zero-argument callable polled between
    sites; returning True raises :class:`MergeAborted`.
    """
    if not sources:
        raise ValueError("no sample sources selected")

    ordered_paths: list[str] = []
    for source in sources:
        if source.file_path not in ordered_paths:
            ordered_paths.append(source.file_path)

    name_by_path = {source.file_path: source.file_name for source in sources}

    # One read handle per backing file.
    cursors: list[_Cursor] = []
    for path in ordered_paths:
        vcf = pysam.VariantFile(path)
        contig_index = {name: idx for idx, name in enumerate(vcf.header.contigs.keys())}
        cursor = _Cursor(
            name=name_by_path[path],
            path=path,
            vcf=vcf,
            record=None,
            contig_index=contig_index,
        )
        cursor.advance()
        cursors.append(cursor)

    # Sample order: first appearance across the ordered source files.
    cohort_samples: list[str] = []
    sample_to_path: dict[str, str] = {}
    for path in ordered_paths:
        for source in sources:
            if source.file_path == path and source.sample not in sample_to_path:
                cohort_samples.append(source.sample)
                sample_to_path[source.sample] = path

    # Output header: rebuild from the first source header (contigs,
    # FORMAT/INFO/FILTER definitions) and replace its sample list with the
    # cohort list. pysam 0.22 has no API to remove samples from a copied
    # header, so replay every meta line on a fresh header instead.
    base_header = pysam.VariantHeader()
    for header_record in cursors[0].vcf.header.records:
        line = str(header_record)
        if line.startswith("#CHROM"):
            continue
        base_header.add_line(line)
    for sample in cohort_samples:
        base_header.add_sample(sample)
    base_header.add_line(
        f'##source=vcf-cohort-merge;samples={len(cohort_samples)};'
        f'source_files={len(ordered_paths)}'
    )

    sites_written = 0
    tmp_path = output_path + ".partial"
    Path(tmp_path).parent.mkdir(parents=True, exist_ok=True)

    writer = pysam.VariantFile(tmp_path, "w", header=base_header)
    try:
        while any(cursor.record is not None for cursor in cursors):
            if is_aborted is not None and is_aborted():
                raise MergeAborted("取用源选择已变更，该在途生成立即失效")

            active = [c for c in cursors if c.record is not None]
            min_key = min(c.sort_key() for c in active)
            chrom = next(
                c.record.chrom for c in active if c.sort_key() == min_key
            )
            pos = min_key[1]

            group = [c for c in active if c.sort_key() == min_key]
            ref = group[0].record.ref

            # Build the ALT union in first-appearance order.
            merged_alts: list[str] = []
            alt_index: dict[str, int] = {}
            for cursor in group:
                record = cursor.record
                for alt in record.alts or ():
                    if alt == ".":
                        continue
                    if alt not in alt_index:
                        alt_index[alt] = len(merged_alts) + 1
                        merged_alts.append(alt)

            output_alleles: tuple[str, ...]
            if merged_alts:
                output_alleles = (ref, *merged_alts)
            else:
                output_alleles = (ref, ".")
            n_out = len(output_alleles)

            out_record = writer.new_record(
                contig=chrom,
                start=pos - 1,
                stop=group[0].record.stop,
                alleles=output_alleles,
            )
            out_record.filter.clear()
            out_record.qual = None
            out_record.id = None

            for sample in cohort_samples:
                owning_path = sample_to_path[sample]
                # The owning cursor is only present when its source file
                # actually carries this coordinate; otherwise the site is
                # missing for this sample and must stay ./. (never 0/0).
                owner = next(
                    (c for c in group if c.path == owning_path), None
                )
                record = owner.record if owner is not None else None

                if record is None or sample not in record.samples:
                    # Site missing from this sample's source file.
                    sample_data = out_record.samples[sample]
                    sample_data["GT"] = (None, None)
                    sample_data["AD"] = tuple([None] * n_out)
                    sample_data["PL"] = tuple(
                        [None] * _triangular_count(n_out)
                    )
                    continue

                src_record = record
                allele_map: list[int | None] = [0]
                for alt in src_record.alts or ():
                    if alt == ".":
                        allele_map.append(None)
                    else:
                        allele_map.append(alt_index[alt])

                sample_data = out_record.samples[sample]
                sample_data["GT"] = _remap_genotype(src_record, sample, allele_map)
                sample_data.phased = src_record.samples[sample].phased
                sample_data["AD"] = _remap_ad(src_record, sample, allele_map, n_out)
                sample_data["PL"] = _remap_pl(src_record, sample, allele_map, n_out)

            writer.write(out_record)
            sites_written += 1

            for cursor in group:
                cursor.advance()
    except BaseException:
        # Abort (including an invalidation noticed mid-merge) or any error:
        # never leave a half-product behind.
        writer.close()
        for cursor in cursors:
            try:
                cursor.vcf.close()
            except Exception:  # noqa: BLE001
                pass
        _safe_unlink(tmp_path)
        raise
    else:
        writer.close()
        for cursor in cursors:
            cursor.vcf.close()

    if is_aborted is not None and is_aborted():
        _safe_unlink(tmp_path)
        raise MergeAborted("取用源选择已变更，该在途生成立即失效")

    os.replace(tmp_path, output_path)
    return {
        "sites": sites_written,
        "samples": len(cohort_samples),
        "output_path": output_path,
    }


def _safe_unlink(path: str) -> None:
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
