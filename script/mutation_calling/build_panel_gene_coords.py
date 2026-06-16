#!/usr/bin/env python3
"""
Build an authoritative gene-coordinate table for the BLYMv2 panel (hg19).

WHY: test.R needs real gene boundaries to map structural-variant breakpoints to
genes. Previously it derived "gene coordinates" from the min/max position of
observed SNVs per gene, which collapses single-mutation genes to a point and
omits genes with no SNVs entirely. This script replaces that with true gene spans
from UCSC refGene, restricted to genes actually targeted by the panel (i.e. whose
span overlaps at least one capture probe in the design BED).

INPUTS:
  --refgene   UCSC refGene.txt(.gz) for hg19. Download once:
              curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
  --bed       BLYMv2 design BED (chr, start, end; 0-based half-open like UCSC).
OUTPUT (TSV, tab-delimited, with header):
  Gene  Chromosome  Gene_Start  Gene_End
  - One row per gene symbol that overlaps the panel.
  - Gene_Start/Gene_End are the union span across all refGene isoforms of that
    symbol on a given chromosome (min txStart .. max txEnd), 1-based start to
    match the coordinate convention used by the rest of the pipeline / VCFs.

USAGE:
  python3 build_panel_gene_coords.py \
      --refgene data/refGene.txt.gz \
      --bed     data/design_BLYMv2_20210721_hg19_sorted.bed \
      --out     data/BLYMv2_gene_coords_hg19.tsv
"""
import argparse
import gzip
import io
import sys
from collections import defaultdict

# Keep only the standard assembly chromosomes (drop *_random / *_alt / chrUn / haps).
STANDARD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, encoding="utf-8")


def load_refgene_gene_spans(refgene_path):
    """Collapse refGene isoforms to one (min txStart, max txEnd) span per (gene, chrom).

    refGene.txt columns (UCSC, with leading bin column):
      0 bin | 1 name(transcript) | 2 chrom | 3 strand | 4 txStart | 5 txEnd
      ... | 12 name2(gene symbol)
    txStart is 0-based; we convert to 1-based for the output start.
    """
    spans = {}  # (gene, chrom) -> [start, end]
    with open_maybe_gz(refgene_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            chrom = f[2]
            if chrom not in STANDARD_CHROMS:
                continue
            try:
                tx_start = int(f[4])
                tx_end = int(f[5])
            except ValueError:
                continue
            gene = f[12].strip()
            if not gene or gene == "":
                continue
            key = (gene, chrom)
            start_1based = tx_start + 1
            if key not in spans:
                spans[key] = [start_1based, tx_end]
            else:
                if start_1based < spans[key][0]:
                    spans[key][0] = start_1based
                if tx_end > spans[key][1]:
                    spans[key][1] = tx_end
    return spans


def load_probes(bed_path):
    """Return {chrom: sorted list of (start_1based, end)} from the design BED (0-based)."""
    probes = defaultdict(list)
    with open_maybe_gz(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            try:
                start = int(f[1]) + 1  # BED 0-based -> 1-based
                end = int(f[2])
            except ValueError:
                continue
            probes[chrom].append((start, end))
    for chrom in probes:
        probes[chrom].sort()
    return probes


def gene_overlaps_panel(chrom, g_start, g_end, probes):
    """True if [g_start, g_end] overlaps any probe interval on chrom."""
    plist = probes.get(chrom)
    if not plist:
        return False
    # Linear scan is fine (a few thousand probes); intervals are sorted by start.
    for p_start, p_end in plist:
        if p_start > g_end:
            break  # sorted: no later probe can overlap
        if p_end >= g_start:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refgene", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    spans = load_refgene_gene_spans(args.refgene)
    probes = load_probes(args.bed)

    rows = []
    for (gene, chrom), (g_start, g_end) in spans.items():
        if gene_overlaps_panel(chrom, g_start, g_end, probes):
            rows.append((gene, chrom, g_start, g_end))

    # Deterministic order: by chrom then start then gene.
    def chrom_key(c):
        c2 = c[3:]  # strip "chr"
        special = {"X": 23, "Y": 24, "M": 25}
        if c2 in special:
            return special[c2]
        try:
            return int(c2)
        except ValueError:
            return 99
    rows.sort(key=lambda r: (chrom_key(r[1]), r[2], r[0]))

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("Gene\tChromosome\tGene_Start\tGene_End\n")
        for gene, chrom, g_start, g_end in rows:
            out.write(f"{gene}\t{chrom}\t{g_start}\t{g_end}\n")

    sys.stderr.write(
        f"Wrote {len(rows)} panel-overlapping genes to {args.out} "
        f"(from {len(spans)} refGene gene/chrom spans).\n"
    )


if __name__ == "__main__":
    main()
