#!/usr/bin/env python3
"""
Create a cCRE BED file for REGENIE gene-based RVAS from ENCODE rE2G predictions.

For each cell type provided, links with Score >= threshold are kept and assigned a
per-cell-type annotation category (cCRE_{label}). The output is a union of all
cell types — one row per (enhancer interval, gene, cell type).

Downstream, each category becomes its own mask in REGENIE (M_cCRE_K562, etc.),
so cell types are tested independently.

Data source — ENCODE portal (no account required):
  wget https://www.encodeproject.org/files/ENCFF950FTI/@@download/ENCFF950FTI.bed.gz  # K562
  wget https://www.encodeproject.org/files/ENCFF055YXP/@@download/ENCFF055YXP.bed.gz  # GM12878
  wget https://www.encodeproject.org/files/ENCFF300OCC/@@download/ENCFF300OCC.bed.gz  # CD34+ (common myeloid progenitor)

File format (header row, tab-separated):
  #chr  start  end  name  class  TargetGene  TargetGeneEnsemblID  ...  Score

Output BED columns (tab-separated, no header):
  chrom   start   end   gene_symbol   category
  chr1    778484  778984  FAM87B      cCRE_K562

Usage:
  uv run python create_ccre_bed.py \\
      --cell-type K562:ENCFF950FTI.bed.gz \\
      --cell-type GM12878:ENCFF055YXP.bed.gz \\
      --cell-type CD34:ENCFF300OCC.bed.gz \\
      --score-threshold 0.5 \\
      --out blood_ccre_gene.bed

After creating the BED, upload to DNAnexus and run ccre_anno.WDL:
  dx upload blood_ccre_gene.bed --destination /resources/
  ./dx_compile.sh
  dx run /workflows/ccre_anno/ccre_anno \\
      -istage-common.step2_chrom_manifest=file-XXX \\
      -istage-common.blood_ccre_gene_bed=file-YYY \\
      --destination /resources/ccre_anno/ -y --brief
"""

import argparse
import sys

import duckdb


def parse_cell_type(value: str) -> tuple[str, str]:
    """Parse 'LABEL:FILE' argument."""
    if ":" not in value:
        raise argparse.ArgumentTypeError(
            f"--cell-type must be LABEL:FILE, got: {value!r}"
        )
    label, path = value.split(":", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError(
            f"--cell-type must be LABEL:FILE with non-empty label and path, got: {value!r}"
        )
    return label, path


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--cell-type",
        action="append",
        metavar="LABEL:FILE",
        required=True,
        dest="cell_types",
        type=parse_cell_type,
        help=(
            "Cell type label and rE2G file path (repeatable). "
            "E.g. --cell-type K562:ENCFF950FTI.bed.gz. "
            "Each label becomes annotation category cCRE_{label}."
        ),
    )
    p.add_argument(
        "--score-threshold",
        type=float,
        default=0.5,
        metavar="T",
        help="Minimum Score required per cell type (default: 0.5).",
    )
    p.add_argument(
        "--tss-bed",
        metavar="FILE",
        help=(
            "3-column TSS BED (no header): chrom, tss_pos, gene_symbol. "
            "When provided, E-G links are filtered by --distance-cap."
        ),
    )
    p.add_argument(
        "--distance-cap",
        type=int,
        default=500_000,
        metavar="BP",
        help=(
            "Max distance (bp) from enhancer midpoint to TSS (default: 500000). "
            "Only applied when --tss-bed is provided."
        ),
    )
    p.add_argument(
        "--out",
        default="blood_ccre_gene.bed",
        help="Output BED file path (default: blood_ccre_gene.bed).",
    )
    return p.parse_args()


def load_cell_type(
    con: duckdb.DuckDBPyConnection,
    path: str,
    label: str,
    threshold: float,
) -> tuple[int, int]:
    """
    Load one rE2G BED file (ENCODE portal format) into table ct_{label}.
    Returns (total_links, links_above_threshold).
    """
    table = f"ct_{label}"
    category = f"cCRE_{label}"
    # ENCODE portal BED files have a header with '#chr' as the first column name.
    # header=true lets DuckDB parse column names; '#chr' must be quoted in SQL.
    con.execute(f"""
        CREATE TABLE {table} AS
        SELECT
            "#chr"                AS chrom,
            CAST("start" AS INT)  AS start,
            CAST("end"   AS INT)  AS end,
            TargetGene            AS gene,
            '{category}'          AS category
        FROM read_csv(
            '{path}',
            delim         = '\t',
            header        = true,
            ignore_errors = true
        )
        WHERE TargetGene IS NOT NULL
          AND TargetGene NOT IN ('', '.')
          AND CAST(Score AS FLOAT) >= {threshold}
    """)
    n_pass = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    # Total before threshold
    con.execute(f"""
        CREATE TEMP TABLE {table}_all AS
        SELECT COUNT(*) AS n
        FROM read_csv(
            '{path}',
            delim         = '\t',
            header        = true,
            ignore_errors = true
        )
        WHERE TargetGene IS NOT NULL AND TargetGene NOT IN ('', '.')
    """)
    n_total = con.execute(f"SELECT n FROM {table}_all").fetchone()[0]
    con.execute(f"DROP TABLE {table}_all")
    return n_total, n_pass


def main():
    args = parse_args()

    if args.score_threshold <= 0 or args.score_threshold >= 1:
        sys.exit("ERROR: --score-threshold must be in (0, 1)")

    labels = [label for label, _ in args.cell_types]
    if len(labels) != len(set(labels)):
        sys.exit("ERROR: duplicate cell type labels")

    con = duckdb.connect()

    for label, path in args.cell_types:
        print(f"Loading {label}: {path}")
        n_total, n_pass = load_cell_type(con, path, label, args.score_threshold)
        n_genes = con.execute(
            f"SELECT COUNT(DISTINCT gene) FROM ct_{label}"
        ).fetchone()[0]
        print(
            f"  {n_total:,} total links | {n_pass:,} with Score >= {args.score_threshold}"
            f" | {n_genes:,} genes"
        )

    # Union all cell types
    union_sql = " UNION ALL ".join(
        f'SELECT chrom, "start", "end", gene, category FROM ct_{label}'
        for label, _ in args.cell_types
    )
    con.execute(f"CREATE TABLE blood_ccre AS {union_sql}")

    n_total = con.execute("SELECT COUNT(*) FROM blood_ccre").fetchone()[0]
    n_genes = con.execute("SELECT COUNT(DISTINCT gene) FROM blood_ccre").fetchone()[0]
    print(f"\nCombined: {n_total:,} E-G links | {n_genes:,} unique genes across all cell types")

    # Per-cell-type summary
    rows = con.execute("""
        SELECT category, COUNT(*) AS n_links, COUNT(DISTINCT gene) AS n_genes
        FROM blood_ccre GROUP BY category ORDER BY category
    """).fetchall()
    for cat, nl, ng in rows:
        print(f"  {cat}: {nl:,} links | {ng:,} genes")

    if args.tss_bed:
        print(
            f"\nApplying distance cap ({args.distance_cap:,} bp) "
            f"using TSS BED: {args.tss_bed}"
        )
        con.execute(f"""
            CREATE TABLE tss AS
            SELECT
                column0              AS chrom,
                CAST(column1 AS INT) AS tss_pos,
                column2              AS gene
            FROM read_csv(
                '{args.tss_bed}',
                delim   = '\t',
                header  = false,
                columns = {{'column0': 'VARCHAR', 'column1': 'INTEGER', 'column2': 'VARCHAR'}}
            )
        """)
        before = con.execute("SELECT COUNT(*) FROM blood_ccre").fetchone()[0]
        con.execute(f"""
            DELETE FROM blood_ccre
            WHERE (chrom, gene) IN (
                SELECT b.chrom, b.gene
                FROM blood_ccre b
                JOIN tss t ON b.chrom = t.chrom AND b.gene = t.gene
                WHERE ABS((b."start" + b."end") / 2 - t.tss_pos) > {args.distance_cap}
            )
        """)
        after = con.execute("SELECT COUNT(*) FROM blood_ccre").fetchone()[0]
        print(f"  Removed {before - after:,} links beyond {args.distance_cap:,} bp from TSS")
        print(f"  {after:,} links remaining")

    # Enhancer width summary
    row = con.execute("""
        SELECT
            ROUND(MEDIAN("end" - "start"), 0) AS median_width,
            MIN("end" - "start")              AS min_width,
            MAX("end" - "start")              AS max_width
        FROM blood_ccre
    """).fetchone()
    median_w, min_w, max_w = row
    print(
        f"\nEnhancer width: median={int(median_w):,} bp  range=[{min_w:,}, {max_w:,}] bp"
    )

    print(f"\nWriting: {args.out}")
    con.execute(f"""
        COPY (
            SELECT chrom, "start", "end", gene, category
            FROM blood_ccre
            ORDER BY chrom, "start", "end", gene, category
        ) TO '{args.out}' (DELIMITER '\t', HEADER false)
    """)
    print(
        f"\nDone. Upload to DNAnexus and run ccre_anno.WDL:\n"
        f"  dx upload {args.out} --destination /resources/\n"
        f"  dx run /workflows/ccre_anno/ccre_anno \\\n"
        f"      -istage-common.step2_chrom_manifest=file-XXX \\\n"
        f"      -istage-common.blood_ccre_gene_bed=file-YYY \\\n"
        f"      --destination /resources/ccre_anno/ -y --brief"
    )


if __name__ == "__main__":
    main()
