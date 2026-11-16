"""Tests for per-file VCF validation rules."""
from __future__ import annotations

import textwrap


def write_vcf(path, body: str, samples_line: str | None = None) -> None:
    header = textwrap.dedent(
        """\
        ##fileformat=VCFv4.2
        ##reference=GRCh38.demo
        ##contig=<ID=chr1,length=100000>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
        """
    )
    if samples_line is None:
        samples_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
    body = textwrap.dedent(body).strip("\n")
    path.write_text(header + samples_line + "\n" + body + "\n")


def test_valid_multiallelic_file(app_env, tmp_path):
    path = tmp_path / "ok.vcf"
    write_vcf(
        path,
        """
        chr1\t100\t.\tA\tG,T\t.\t.\t.\tGT:AD:PL\t0/1:10,5,2:0,9,200,12,210,215
        chr1\t200\t.\tC\tG\t.\t.\t.\tGT:AD:PL\t1/1:1,9:220,12,0
        """,
    )
    result = app_env["validation"].validate_vcf(str(path))
    assert result.valid, result.error
    assert result.samples == ["S1"]
    assert result.record_count == 2
    assert result.reference_version == "GRCh38.demo"


def test_unsorted_file_rejected(app_env, tmp_path):
    path = tmp_path / "unsorted.vcf"
    write_vcf(
        path,
        """
        chr1\t200\t.\tC\tG\t.\t.\t.\tGT\t0/0
        chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1
        """,
    )
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "坐标" in result.error


def test_duplicate_coordinate_rejected(app_env, tmp_path):
    path = tmp_path / "dup.vcf"
    write_vcf(
        path,
        """
        chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1
        chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1
        """,
    )
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "同坐标仅一条记录" in result.error


def test_non_acgt_allele_rejected(app_env, tmp_path):
    path = tmp_path / "allele.vcf"
    write_vcf(path, "chr1\t100\t.\tA\tGN\t.\t.\t.\tGT\t0/1\n")
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "ACGT" in result.error


def test_forbidden_format_field_rejected(app_env, tmp_path):
    path = tmp_path / "fmt.vcf"
    write_vcf(
        path,
        "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:15\n",
    )
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "FORMAT" in result.error


def test_missing_gt_field_rejected(app_env, tmp_path):
    path = tmp_path / "nogt.vcf"
    write_vcf(path, "chr1\t100\t.\tA\tG\t.\t.\t.\tAD\t10,5\n")
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "GT" in result.error


def test_non_diploid_rejected(app_env, tmp_path):
    path = tmp_path / "triploid.vcf"
    write_vcf(path, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1/1\n")
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "二倍体" in result.error


def test_missing_diploid_genotype_allowed(app_env, tmp_path):
    path = tmp_path / "missing.vcf"
    write_vcf(path, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t./.:.:.\n")
    result = app_env["validation"].validate_vcf(str(path))
    assert result.valid, result.error


def test_empty_file_rejected(app_env, tmp_path):
    path = tmp_path / "empty.vcf"
    write_vcf(path, "")
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid


def test_undeclared_contig_rejected(app_env, tmp_path):
    path = tmp_path / "contig.vcf"
    write_vcf(path, "chr9\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")
    result = app_env["validation"].validate_vcf(str(path))
    assert not result.valid
    assert "contig" in result.error
