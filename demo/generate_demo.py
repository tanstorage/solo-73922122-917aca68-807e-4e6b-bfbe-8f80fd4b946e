#!/usr/bin/env python3
"""Generate deterministic multi-allelic diploid demo VCF files.

Run with the project virtualenv, e.g.::

    python demo/generate_demo.py demo

Three valid, coordinate-sorted VCFs sharing one reference version are written
plus one deliberately invalid VCF used to demonstrate whole-file rejection:

* batch_a.vcf : samples S1, S2 (S2 is duplicated in batch_b)
* batch_b.vcf : samples S2, S3
* batch_c.vcf : samples S3, S4 (S3 duplicated in batch_b)
* bad_ref_conflict.vcf : same coordinate as batch_a but a different REF

The trio covers:
  - multiallelic sites with 2-3 ALT alleles,
  - ALT unions whose order differs across files,
  - missing sites in a file (must stay ./. in the cohort, never 0/0),
  - no-ALT (monomorphic) sites,
  - phased and unphased genotypes,
  - fully missing genotypes and partially missing AD/PL.
"""
from __future__ import annotations

import sys
from pathlib import Path

HEADER_TEMPLATE = """##fileformat=VCFv4.2
##reference=GRCh38.demo.1.0
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""


def write_vcf(path: Path, samples: list[str], rows: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(HEADER_TEMPLATE.format(samples="\t".join(samples)))
        handle.write("\n".join(rows))
        handle.write("\n")


def main(out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- batch_a: S1, S2 --------------------------------------------------
    rows_a = [
        # Site 1: A -> G,T in file A
        "chr1\t100\t.\tA\tG,T\t.\t.\t.\tGT:AD:PL\t"
        "0/1:12,8,0:0,9,210,14,215,215\t"
        "1/2:3,10,9:240,220,12,230,18,0",
        # Site 2: monomorphic (no ALT)
        "chr1\t250\t.\tC\t.\t.\t.\t.\tGT:AD:PL\t"
        "0/0:20:0,30,250\t"
        "0/0:18:0,28,240",
        # Site 3: C -> G only here; S2 fully missing
        "chr1\t400\t.\tC\tG\t.\t.\t.\tGT:AD:PL\t"
        "1/1:2,15:220,25,0\t"
        "./.:.:.",
        # Site 4 on chr2: T -> A,C (phased genotype)
        "chr2\t100\t.\tT\tA,C\t.\t.\t.\tGT:AD:PL\t"
        "0|1:11,9,1:0,10,200,13,210,220\t"
        "0/0:19,0,0:0,30,250,32,255,255",
    ]
    write_vcf(out / "batch_a.vcf", ["S1", "S2"], rows_a)

    # ---- batch_b: S2, S3 --------------------------------------------------
    rows_b = [
        # Site 1: same REF, ALT order T then G (union order must follow first
        # appearance in file order: G,T from batch_a)
        "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:AD:PL\t"
        "0/1:10,7,2:0,12,180,16,200,210\t"
        "2/2:1,1,17:250,240,230,220,210,0",
        # Site 3 absent -> S2/S3 must be ./. at chr1:400 in cohort
        # Site 4: introduces a new ALT G -> union A,C,G
        "chr2\t100\t.\tT\tG\t.\t.\t.\tGT:AD:PL\t"
        "0/1:9,8:0,11,190\t"
        "1/1:0,14:230,20,0",
        # New site only in B
        "chr2\t300\t.\tG\tA,T\t.\t.\t.\tGT:AD:PL\t"
        "0/2:8,1,7:0,160,14,200,210,180\t"
        "./.:.:.",
    ]
    write_vcf(out / "batch_b.vcf", ["S2", "S3"], rows_b)

    # ---- batch_c: S3, S4 --------------------------------------------------
    rows_c = [
        # Site 1: A -> G only; contributes no new ALT to the union
        "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:PL\t"
        "0/0:16,1:0,12,190\t"
        "0/1:10,6:0,15,200",
        # Site 3 present with same REF as batch_a
        "chr1\t400\t.\tC\tG\t.\t.\t.\tGT:AD:PL\t"
        "0/1:13,7:0,8,180\t"
        "1/1:1,17:225,20,0",
        # Site 4: yet another ALT A already in union; adds nothing new
        "chr2\t100\t.\tT\tA\t.\t.\t.\tGT:AD:PL\t"
        "1/1:2,16:230,18,0\t"
        "./.:.:.",
        "chr2\t300\t.\tG\tA\t.\t.\t.\tGT:AD:PL\t"
        "0/1:12,9:0,10,200\t"
        "1/1:0,20:240,22,0",
    ]
    write_vcf(out / "batch_c.vcf", ["S3", "S4"], rows_c)

    # ---- deliberately invalid file: REF conflict at chr1:100 -------------
    bad_rows = [
        "chr1\t100\t.\tAA\tG\t.\t.\t.\tGT:AD:PL\t"
        "0/1:10,5:0,9,200\t"
        "0/0:14:0,20,230",
    ]
    write_vcf(out / "bad_ref_conflict.vcf", ["S5", "S6"], bad_rows)

    print(f"demo files written to {out.resolve()}")
    for name in ("batch_a.vcf", "batch_b.vcf", "batch_c.vcf", "bad_ref_conflict.vcf"):
        print(" -", name)


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "demo"
    main(target)
