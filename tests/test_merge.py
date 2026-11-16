"""Tests for the cohort merge engine: ALT union and GT/AD/PL remapping."""
from __future__ import annotations

import textwrap

import pysam
import pytest

from app.merge import SampleSource, merge_vcfs, _pl_index, MergeAborted


HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##reference=GRCh38.demo
    ##contig=<ID=chr1,length=100000>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
    """
)


def vcf_with(path, samples: list[str], rows: list[str]) -> None:
    header = HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    path.write_text(header + "\n".join(rows) + "\n")


def read_output(path):
    vcf = pysam.VariantFile(str(path))
    records = list(vcf)
    samples = list(vcf.header.samples)
    vcf.close()
    return samples, records


def test_alt_union_and_gt_ad_pl_remap(tmp_path):
    f1 = tmp_path / "a.vcf"
    f2 = tmp_path / "b.vcf"
    # File A: A -> G,T ; File B: A -> T,G (different ALT order)
    vcf_with(
        f1, ["S1"],
        ["chr1\t100\t.\tA\tG,T\t.\t.\t.\tGT:AD:PL\t0/1:10,5,2:0,9,200,12,210,215"],
    )
    vcf_with(
        f2, ["S2"],
        ["chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:AD:PL\t0/2:8,1,6:0,180,11,200,210,205"],
    )
    out = tmp_path / "cohort.vcf"
    merge_vcfs(
        [
            SampleSource("S1", str(f1), "a.vcf"),
            SampleSource("S2", str(f2), "b.vcf"),
        ],
        str(out),
    )
    samples, records = read_output(out)
    assert samples == ["S1", "S2"]
    assert len(records) == 1
    record = records[0]
    # Union order follows first appearance across ordered files: G then T.
    assert record.alleles == ("A", "G", "T")

    s1 = record.samples["S1"]
    assert s1["GT"] == (0, 1)
    assert tuple(s1["AD"]) == (10, 5, 2)
    assert tuple(s1["PL"]) == (0, 9, 200, 12, 210, 215)

    # File B used T=idx1, G=idx2; in the union G=1, T=2, so 0/2 -> 0/1.
    s2 = record.samples["S2"]
    assert s2["GT"] == (0, 1)
    # B had no value for union allele T slot coming from B's own idx1 =>
    # AD remap: B AD=(8,1,6) for (REF, T, G) -> union (REF=8, G=6, T=1)
    assert tuple(s2["AD"]) == (8, 6, 1)
    # PL remap for diploid triangular order. Source B PL layout for alleles
    # (REF,T,G): (00)=0,(01)=180,(11)=11,(02)=200,(12)=210,(22)=205.
    # Union order is (REF,G,T): B's (REF,G)=source 02 -> union 01.
    pl = tuple(s2["PL"])
    assert len(pl) == 6
    assert pl[_pl_index(0, 1, 3)] == 200
    # B's homozygous T (source 11, value 11) maps to union 22.
    assert pl[_pl_index(2, 2, 3)] == 11
    assert pl[0] == 0


def test_new_allele_keeps_missing_values(tmp_path):
    f1 = tmp_path / "a.vcf"
    f2 = tmp_path / "b.vcf"
    vcf_with(
        f1, ["S1"],
        ["chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t0/1:10,5:0,9,200"],
    )
    vcf_with(
        f2, ["S2"],
        ["chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD:PL\t0/1:8,4:0,7,190"],
    )
    out = tmp_path / "cohort.vcf"
    merge_vcfs(
        [
            SampleSource("S1", str(f1), "a.vcf"),
            SampleSource("S2", str(f2), "b.vcf"),
        ],
        str(out),
    )
    _, records = read_output(out)
    record = records[0]
    assert record.alleles == ("A", "G", "T")
    s1 = record.samples["S1"]
    # T was introduced by file B: S1's T-related values must be missing.
    assert tuple(s1["AD"]) == (10, 5, None)
    pl = tuple(s1["PL"])
    # Source PL layout (REF,G): 00=0, 01=9, 11=200 -> union indices 00, 01, 11.
    assert pl[0] == 0 and pl[1] == 9 and pl[_pl_index(1, 1, 3)] == 200
    # Union positions involving the new T allele (02, 12, 22) stay missing.
    assert pl[_pl_index(0, 2, 3)] is None
    assert pl[_pl_index(1, 2, 3)] is None
    assert pl[_pl_index(2, 2, 3)] is None
    s2 = record.samples["S2"]
    assert tuple(s2["AD"]) == (8, None, 4)


def test_missing_site_is_missing_not_homref(tmp_path):
    f1 = tmp_path / "a.vcf"
    f2 = tmp_path / "b.vcf"
    vcf_with(
        f1, ["S1"],
        [
            "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t0/1:10,5:0,9,200",
            "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:AD:PL\t1/1:1,9:220,12,0",
        ],
    )
    vcf_with(
        f2, ["S2"],
        ["chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t0/0:14,0:0,20,230"],
    )
    out = tmp_path / "cohort.vcf"
    merge_vcfs(
        [
            SampleSource("S1", str(f1), "a.vcf"),
            SampleSource("S2", str(f2), "b.vcf"),
        ],
        str(out),
    )
    samples, records = read_output(out)
    by_pos = {r.pos: r for r in records}
    assert set(by_pos) == {100, 200}
    missing_sample = by_pos[200].samples["S2"]
    assert missing_sample["GT"] == (None, None)
    assert all(v is None for v in missing_sample["AD"])
    assert all(v is None for v in missing_sample["PL"])
    # And it must not have been turned into 0/0.
    text = out.read_text()
    line200 = [line for line in text.splitlines() if line.startswith("chr1\t200")][0]
    columns = line200.split("\t")
    assert columns[9 + samples.index("S2")].startswith("./.")
    # The present sample keeps its real genotype.
    assert columns[9 + samples.index("S1")].startswith("1/1")


def test_duplicate_sample_uses_one_original_file_only(tmp_path):
    f1 = tmp_path / "a.vcf"
    f2 = tmp_path / "b.vcf"
    vcf_with(
        f1, ["DUP"],
        ["chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t0/1:10,5:0,9,200"],
    )
    vcf_with(
        f2, ["DUP"],
        [
            "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t1/1:1,9:220,12,0",
            "chr1\t300\t.\tT\tA\t.\t.\t.\tGT:AD:PL\t0/1:7,4:0,6,180",
        ],
    )
    out_a = tmp_path / "from_a.vcf"
    merge_vcfs([SampleSource("DUP", str(f1), "a.vcf")], str(out_a))
    _, records_a = read_output(out_a)
    assert len(records_a) == 1  # no site from B may be spliced in
    assert records_a[0].samples["DUP"]["GT"] == (0, 1)

    out_b = tmp_path / "from_b.vcf"
    merge_vcfs([SampleSource("DUP", str(f2), "b.vcf")], str(out_b))
    _, records_b = read_output(out_b)
    assert len(records_b) == 2
    assert records_b[0].samples["DUP"]["GT"] == (1, 1)


def test_phasing_is_preserved(tmp_path):
    f1 = tmp_path / "a.vcf"
    vcf_with(
        f1, ["S1"],
        ["chr1\t100\t.\tA\tG,T\t.\t.\t.\tGT:AD:PL\t1|2:3,10,9:240,220,12,230,18,0"],
    )
    out = tmp_path / "cohort.vcf"
    merge_vcfs([SampleSource("S1", str(f1), "a.vcf")], str(out))
    _, records = read_output(out)
    sample = records[0].samples["S1"]
    assert sample.phased is True
    assert sample["GT"] == (1, 2)


def test_no_alt_site_written_with_dot_alt(tmp_path):
    f1 = tmp_path / "a.vcf"
    f2 = tmp_path / "b.vcf"
    vcf_with(f1, ["S1"], ["chr1\t100\t.\tC\t.\t.\t.\t.\tGT\t0/0"])
    vcf_with(
        f2, ["S2"],
        ["chr1\t100\t.\tC\tG\t.\t.\t.\tGT:AD:PL\t0/1:8,4:0,7,190"],
    )
    out = tmp_path / "cohort.vcf"
    merge_vcfs(
        [
            SampleSource("S1", str(f1), "a.vcf"),
            SampleSource("S2", str(f2), "b.vcf"),
        ],
        str(out),
    )
    _, records = read_output(out)
    record = records[0]
    # Union contains the G allele contributed by file B.
    assert record.alleles == ("C", "G")
    assert record.samples["S1"]["GT"] == (0, 0)
    assert record.samples["S2"]["GT"] == (0, 1)
    assert "C\tG" in out.read_text()


def test_merge_aborted_removes_partial_output(tmp_path):
    f1 = tmp_path / "a.vcf"
    vcf_with(
        f1, ["S1"],
        [f"chr1\t{100 + i}\t.\tA\tG\t.\t.\t.\tGT\t0/1" for i in range(3)],
    )
    out = tmp_path / "cohort.vcf"
    with pytest.raises(MergeAborted):
        merge_vcfs(
            [SampleSource("S1", str(f1), "a.vcf")],
            str(out),
            is_aborted=lambda: True,
        )
    assert not out.exists()
    assert not (tmp_path / "cohort.vcf.partial").exists()


def test_pl_index_layout():
    # 3 alleles -> 6 genotypes in VCF triangular order.
    order = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]
    for idx, (i, j) in enumerate(order):
        assert _pl_index(i, j, 3) == idx
