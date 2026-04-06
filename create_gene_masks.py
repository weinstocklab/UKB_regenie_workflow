#!/usr/bin/env python3
"""
Create REGENIE gene-based test annotation files from AlphaMissense hg38 predictions
at multiple am_pathogenicity score thresholds.

Produces three files in --out-dir:
  anno_file.tsv  -- variant_id  gene  category
  set_list.tsv   -- gene  chrom  start_pos  variant1,variant2,...
  mask_def.tsv   -- mask_name  annotation_categories

AlphaMissense hg38 download:
  https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz

Variant IDs are constructed as PREFIX:CHROM:POS:REF:ALT (e.g. DRAGEN:chr1:12345:A:G).
Use --variant-id-prefix to match the IDs in your WGS pvar files (default: DRAGEN).

Gene assignment (choose one, in priority order):
  --uniprot-mapping  UniProt accession -> gene symbol TSV (recommended).
                     Download from UniProt REST API:
                       curl -s "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606&fields=accession,gene_names&format=tsv&compressed=true" \\
                           -o uniprot_human_gene_names.tsv.gz
                     Format: tab-separated with header (Entry, Gene Names).
                     "Gene Names" is space-separated; the first token is the primary HGNC symbol.
  --gene-mapping     TSV with columns transcript_id, gene_name (e.g. from GENCODE/Ensembl).
                     Transcript versions are stripped before matching (ENST0001.4 -> ENST0001).
                     Generate with:
                       zcat gencode.v45.annotation.gtf.gz \\
                         | awk '$3=="transcript"' \\
                         | grep -oP 'transcript_id "[^"]+"|gene_name "[^"]+"' \\
                         | paste - - \\
                         | sed 's/transcript_id "//;s/gene_name "//;s/"//g' \\
                         > transcript_gene.tsv
  (default)          Use uniprot_id from AlphaMissense as the gene identifier.
                     Works without external files but uses UniProt accessions, not gene symbols.

Multiple thresholds:
  Each threshold creates an annotation category AM_<threshold> (e.g. AM_0p34 for 0.34).
  Each (variant, gene) pair is assigned exactly ONE category — the strictest threshold
  it qualifies for (REGENIE requires at most one annotation per variant-gene pair).
  A variant with score 0.7 and thresholds [0.34, 0.564] gets one anno_file row:
    chr1-xxx  GENE  AM_0p564      ← strictest passing threshold
  A variant with score 0.45 gets:
    chr1-xxx  GENE  AM_0p34       ← only threshold it passes
  Masks include all categories at and above the threshold:
    M_AM_0p34     AM_0p34,AM_0p564    ← captures score >= 0.34
    M_AM_0p564    AM_0p564            ← captures score >= 0.564

  Known AlphaMissense category boundaries (from this file):
    likely_benign     : AM < 0.34
    ambiguous (VUS)   : 0.34 <= AM < 0.564
    likely_pathogenic : AM >= 0.564

Optional pLoF input:
  --plof  Tab-separated file with columns variant_id, gene_name (no header).
          Variant IDs must use the same PREFIX:CHROM:POS:REF:ALT format.

Usage:
  uv run python create_gene_masks.py \\
      --alphamissense AlphaMissense_hg38.tsv.gz \\
      --uniprot-mapping uniprot_human_gene_names.tsv.gz \\
      --thresholds 0.34 0.564 \\
      --variant-id-prefix DRAGEN \\
      --out-dir gene_masks/
"""

import argparse
import re
import sys

import duckdb


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--alphamissense",
        required=True,
        help="AlphaMissense hg38 TSV or TSV.gz.",
    )
    p.add_argument(
        "--uniprot-mapping",
        metavar="FILE",
        help=(
            "UniProt TSV (with header: Entry, Gene Names) mapping accessions to gene symbols. "
            "Download: curl -s 'https://rest.uniprot.org/uniprotkb/stream?"
            "query=organism_id:9606&fields=accession,gene_names&format=tsv&compressed=true' "
            "-o uniprot_human_gene_names.tsv.gz"
        ),
    )
    p.add_argument(
        "--gene-mapping",
        metavar="FILE",
        help=(
            "Two-column TSV (no header): transcript_id, gene_name. "
            "Used instead of --uniprot-mapping when transcript-level mapping is preferred."
        ),
    )
    p.add_argument(
        "--plof",
        metavar="FILE",
        help=(
            "Two-column TSV (no header): variant_id (CHROM-POS-REF-ALT), gene_name. "
            "Adds pLoF annotation category and pLoF-containing masks."
        ),
    )
    p.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[0.34, 0.564],
        metavar="T",
        help=(
            "AM score thresholds (inclusive lower bound per category). "
            "Default: 0.34 0.564  (benign/VUS and VUS/pathogenic boundaries)."
        ),
    )
    p.add_argument(
        "--variant-id-prefix",
        default="DRAGEN",
        metavar="PREFIX",
        help=(
            "Prefix prepended to variant IDs: PREFIX:CHROM:POS:REF:ALT. "
            "Must match the IDs in your WGS pvar files (default: DRAGEN)."
        ),
    )
    p.add_argument(
        "--out-dir",
        default=".",
        help="Output directory (default: current directory).",
    )
    return p.parse_args()


def threshold_to_name(t: float) -> str:
    """Convert a float threshold to a safe category name suffix, e.g. 0.564 -> '0p564'."""
    s = f"{t:g}"          # '0.564' or '0.34'
    s = s.replace(".", "p")  # '0p564'
    return s


def load_am_with_gene(
    con,
    am_path: str,
    uniprot_mapping_path: str | None,
    gene_mapping_path: str | None,
    variant_id_prefix: str = "DRAGEN",
) -> int:
    """
    Load AlphaMissense, map to gene symbols, take MAX score per (variant, gene).
    Priority: uniprot_mapping > gene_mapping (transcript-based) > uniprot_id as-is.
    Creates table `am_by_gene` with columns: variant_id, gene, max_score.
    Returns row count.
    """
    # The AlphaMissense file has leading '#...' comment lines then a '#CHROM' header.
    # Use comment='#' + explicit names to skip all comment/header lines.
    print(f"Loading AlphaMissense: {am_path}")
    print(f"  Variant ID format: {variant_id_prefix}:CHROM:POS:REF:ALT")
    con.execute(f"""
        CREATE OR REPLACE TABLE am_raw AS
        SELECT
            '{variant_id_prefix}' || ':' || CHROM || ':' || CAST(POS AS VARCHAR) || ':' || REF || ':' || ALT
                AS variant_id,
            -- strip transcript version suffix (ENST00001.4 -> ENST00001)
            regexp_replace(transcript_id, '\\.\\d+$', '') AS transcript_id,
            uniprot_id,
            CAST(am_pathogenicity AS FLOAT) AS am_pathogenicity
        FROM read_csv(
            '{am_path}',
            delim = '\t',
            header = false,
            comment = '#',
            names = ['CHROM','POS','REF','ALT','genome','uniprot_id',
                     'transcript_id','protein_variant','am_pathogenicity','am_class']
        )
    """)
    n_raw = con.execute("SELECT COUNT(*) FROM am_raw").fetchone()[0]
    print(f"  {n_raw:,} total (variant, transcript) rows")

    if uniprot_mapping_path:
        # UniProt REST API format: header row "Entry\tGene Names"
        # "Gene Names" is space-separated; first token is the primary HGNC symbol.
        print(f"  Mapping uniprot_id -> gene symbol via: {uniprot_mapping_path}")
        con.execute(f"""
            CREATE OR REPLACE TABLE uniprot_map AS
            SELECT
                "Entry" AS uniprot_id,
                split_part("Gene Names", ' ', 1) AS gene_name
            FROM read_csv(
                '{uniprot_mapping_path}',
                delim = '\t',
                header = true
            )
            WHERE "Gene Names" IS NOT NULL AND "Gene Names" != ''
        """)
        n_map = con.execute("SELECT COUNT(*) FROM uniprot_map").fetchone()[0]
        print(f"  Loaded {n_map:,} UniProt -> gene symbol entries")

        con.execute("""
            CREATE OR REPLACE TABLE am_by_gene AS
            SELECT
                a.variant_id,
                u.gene_name AS gene,
                MAX(a.am_pathogenicity) AS max_score
            FROM am_raw a
            JOIN uniprot_map u ON a.uniprot_id = u.uniprot_id
            GROUP BY a.variant_id, u.gene_name
        """)
        n_unmatched = con.execute("""
            SELECT COUNT(DISTINCT uniprot_id) FROM am_raw
            WHERE uniprot_id NOT IN (SELECT uniprot_id FROM uniprot_map)
        """).fetchone()[0]
        if n_unmatched > 0:
            print(f"  Warning: {n_unmatched:,} AM uniprot_ids had no gene mapping")

    elif gene_mapping_path:
        print(f"  Mapping transcripts -> genes via: {gene_mapping_path}")
        con.execute(f"""
            CREATE OR REPLACE TABLE gene_map AS
            SELECT
                regexp_replace(column0, '\\.\\d+$', '') AS transcript_id,
                column1 AS gene_name
            FROM read_csv(
                '{gene_mapping_path}',
                delim = '\t',
                header = false,
                columns = {{'column0': 'VARCHAR', 'column1': 'VARCHAR'}}
            )
        """)
        n_mapped = con.execute("SELECT COUNT(DISTINCT transcript_id) FROM gene_map").fetchone()[0]
        print(f"  Gene mapping covers {n_mapped:,} unique transcripts")

        con.execute("""
            CREATE OR REPLACE TABLE am_by_gene AS
            SELECT
                a.variant_id,
                g.gene_name AS gene,
                MAX(a.am_pathogenicity) AS max_score
            FROM am_raw a
            JOIN gene_map g ON a.transcript_id = g.transcript_id
            GROUP BY a.variant_id, g.gene_name
        """)
        n_unmatched = con.execute("""
            SELECT COUNT(DISTINCT transcript_id) FROM am_raw
            WHERE transcript_id NOT IN (SELECT transcript_id FROM gene_map)
        """).fetchone()[0]
        if n_unmatched > 0:
            print(f"  Warning: {n_unmatched:,} AM transcripts had no gene mapping")

    else:
        print("  No gene mapping provided — using uniprot_id as gene identifier")
        con.execute("""
            CREATE OR REPLACE TABLE am_by_gene AS
            SELECT
                variant_id,
                uniprot_id AS gene,
                MAX(am_pathogenicity) AS max_score
            FROM am_raw
            GROUP BY variant_id, uniprot_id
        """)

    n = con.execute("SELECT COUNT(*) FROM am_by_gene").fetchone()[0]
    n_genes = con.execute("SELECT COUNT(DISTINCT gene) FROM am_by_gene").fetchone()[0]
    n_variants = con.execute("SELECT COUNT(DISTINCT variant_id) FROM am_by_gene").fetchone()[0]
    print(f"  {n:,} (variant, gene) pairs | {n_variants:,} unique variants | {n_genes:,} genes")
    return n


def build_am_annotations(con, thresholds: list[float]) -> int:
    """
    For each threshold T, annotate qualifying (variant, gene) pairs with category AM_<T>.
    A variant with score 0.7 and thresholds [0.34, 0.564] gets two rows.
    Creates table `am_anno` with columns: variant_id, gene, category.
    Returns total row count.
    """
    print(f"\nBuilding AM annotations at thresholds: {thresholds}")
    min_thresh = min(thresholds)

    # Build a VALUES table for the thresholds
    values_rows = ", ".join(
        f"({t}, 'AM_{threshold_to_name(t)}')" for t in thresholds
    )
    con.execute(f"""
        CREATE OR REPLACE TABLE thresh_table AS
        SELECT * FROM (VALUES {values_rows}) t(threshold, category)
    """)

    # Assign each (variant, gene) exactly ONE category: the strictest threshold passed.
    # REGENIE requires at most one annotation per (variant, gene) pair.
    # Masks then aggregate categories: M_AM_0p34 = AM_0p34,AM_0p564 to capture both bins.
    con.execute(f"""
        CREATE OR REPLACE TABLE am_anno AS
        SELECT a.variant_id, a.gene, t.category
        FROM am_by_gene a
        CROSS JOIN thresh_table t
        WHERE a.max_score >= t.threshold
          AND a.max_score >= {min_thresh}
        QUALIFY t.threshold = MAX(t.threshold) OVER (PARTITION BY a.variant_id, a.gene)
    """)
    n = con.execute("SELECT COUNT(*) FROM am_anno").fetchone()[0]

    # Breakdown by category
    rows = con.execute("""
        SELECT category, COUNT(*) AS n_pairs, COUNT(DISTINCT variant_id) AS n_variants
        FROM am_anno
        GROUP BY category
        ORDER BY category
    """).fetchall()
    for cat, n_pairs, n_var in rows:
        print(f"  {cat}: {n_pairs:,} variant-gene pairs  ({n_var:,} unique variants)")
    return n


def load_plof(con, plof_path: str) -> int:
    """Load pLoF file (2-col TSV: variant_id, gene_name, no header). Returns row count."""
    print(f"\nLoading pLoF annotations: {plof_path}")
    con.execute(f"""
        CREATE OR REPLACE TABLE plof AS
        SELECT DISTINCT
            column0 AS variant_id,
            column1 AS gene,
            'pLoF'  AS category
        FROM read_csv(
            '{plof_path}',
            delim = '\t',
            header = false,
            columns = {{'column0': 'VARCHAR', 'column1': 'VARCHAR'}}
        )
    """)
    n = con.execute("SELECT COUNT(*) FROM plof").fetchone()[0]
    print(f"  {n:,} pLoF variant-gene pairs")
    return n


def build_combined(con, has_plof: bool) -> int:
    """
    Merge pLoF (if present) and AM annotations.
    pLoF takes precedence: an AM annotation for the same (variant, gene) is dropped.
    Creates table `anno` with columns: variant_id, gene, category.
    """
    print("\nBuilding combined annotation table...")
    if has_plof:
        con.execute("""
            CREATE OR REPLACE TABLE anno AS
            SELECT variant_id, gene, category FROM plof
            UNION ALL
            SELECT a.variant_id, a.gene, a.category
            FROM am_anno a
            WHERE NOT EXISTS (
                SELECT 1 FROM plof p
                WHERE p.variant_id = a.variant_id AND p.gene = a.gene
            )
        """)
    else:
        con.execute("""
            CREATE OR REPLACE TABLE anno AS
            SELECT variant_id, gene, category FROM am_anno
        """)

    n = con.execute("SELECT COUNT(*) FROM anno").fetchone()[0]
    n_genes = con.execute("SELECT COUNT(DISTINCT gene) FROM anno").fetchone()[0]
    cats = con.execute("""
        SELECT category, COUNT(*) AS n FROM anno GROUP BY category ORDER BY category
    """).fetchall()
    print(f"  {n:,} total rows | {n_genes:,} genes")
    for cat, cnt in cats:
        print(f"    {cat}: {cnt:,}")
    return n


def write_anno_file(con, out_dir: str):
    path = f"{out_dir}/anno_file.tsv"
    n = con.execute("SELECT COUNT(*) FROM anno").fetchone()[0]
    con.execute(f"""
        COPY (
            SELECT DISTINCT variant_id, gene, category
            FROM anno
            ORDER BY gene, variant_id, category
        ) TO '{path}' (DELIMITER '\t', HEADER false)
    """)
    print(f"Written: {path}  ({n:,} rows)")


def write_set_list(con, out_dir: str):
    # REGENIE format: setID  chrom  start_pos  variant_list
    # Variant ID format: PREFIX:CHROM:POS:REF:ALT  (e.g. DRAGEN:chr9:10001:T:A)
    # split_part is 1-based: part 2 = CHROM, part 3 = POS
    con.execute("""
        CREATE OR REPLACE TABLE set_list_tbl AS
        SELECT
            gene,
            regexp_replace(split_part(variant_id, ':', 2), '^chr', '') AS chrom,
            MIN(CAST(split_part(variant_id, ':', 3) AS INTEGER)) AS start_pos,
            string_agg(DISTINCT variant_id, ',') AS variant_list
        FROM anno
        GROUP BY gene, regexp_replace(split_part(variant_id, ':', 2), '^chr', '')
    """)
    n = con.execute("SELECT COUNT(*) FROM set_list_tbl").fetchone()[0]
    path = f"{out_dir}/set_list.tsv"
    con.execute(f"""
        COPY (
            SELECT gene, chrom, start_pos, variant_list
            FROM set_list_tbl
            WHERE TRY_CAST(chrom AS INTEGER) BETWEEN 1 AND 22
            ORDER BY CAST(chrom AS INTEGER), start_pos
        ) TO '{path}' (DELIMITER '\t', HEADER false)
    """)
    n_written = con.execute(
        "SELECT COUNT(*) FROM set_list_tbl WHERE TRY_CAST(chrom AS INTEGER) BETWEEN 1 AND 22"
    ).fetchone()[0]
    n_skipped = n - n_written
    if n_skipped > 0:
        print(f"  Skipped {n_skipped:,} non-autosomal gene/chrom entries (X, Y, M, etc.)")
    print(f"Written: {path}  ({n_written:,} autosomal genes/sets)")


def write_mask_def(thresholds: list[float], has_plof: bool, out_dir: str):
    """
    Auto-generate mask definitions.
    Thresholds are sorted ascending so masks go from broader to stricter.
    """
    thresholds_asc = sorted(thresholds)
    cats = [f"AM_{threshold_to_name(t)}" for t in thresholds_asc]

    # For each threshold, the mask includes that category AND all stricter (higher) ones,
    # because stricter-threshold variants are binned into the stricter category only.
    # e.g. M_AM_0p34 = AM_0p34,AM_0p564 captures all variants with score >= 0.34.
    masks = []
    if has_plof:
        masks.append(("M_pLoF", "pLoF"))
        for i, cat in enumerate(cats):
            included = "pLoF," + ",".join(cats[i:])
            masks.append((f"M_pLoF_{cat}", included))
    for i, cat in enumerate(cats):
        included = ",".join(cats[i:])
        masks.append((f"M_{cat}", included))

    path = f"{out_dir}/mask_def.tsv"
    with open(path, "w") as f:
        for name, definition in masks:
            f.write(f"{name}\t{definition}\n")
    print(f"Written: {path}")
    print("\nMask definitions:")
    for name, definition in masks:
        print(f"  {name:<30}  {definition}")


def main():
    args = parse_args()

    # Validate thresholds
    if any(t <= 0 or t >= 1 for t in args.thresholds):
        sys.exit("ERROR: all thresholds must be in (0, 1)")
    thresholds = sorted(set(args.thresholds))

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    con = duckdb.connect()

    load_am_with_gene(con, args.alphamissense, args.uniprot_mapping, args.gene_mapping, args.variant_id_prefix)
    build_am_annotations(con, thresholds)

    has_plof = args.plof is not None
    if has_plof:
        load_plof(con, args.plof)

    build_combined(con, has_plof)

    print()
    write_anno_file(con, args.out_dir)
    write_set_list(con, args.out_dir)
    write_mask_def(thresholds, has_plof, args.out_dir)

    print("\nDone. Upload to DNAnexus and pass file IDs to regenie_gene_bt/qt.WDL:")
    print("  --anno-file   anno_file.tsv")
    print("  --set-list    set_list.tsv")
    print("  --mask-def    mask_def.tsv")


if __name__ == "__main__":
    main()
