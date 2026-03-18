version 1.0

# ccre_anno.WDL
#
# On-platform workflow to intersect ENCODE cCRE enhancer-gene intervals with
# UKB WGS variant positions (pvar files on DNAnexus), producing the three
# annotation files required by regenie_gene_bt/qt.WDL:
#   anno_file.tsv  -- variant_id  gene  category
#   set_list.tsv   -- gene  chrom  start_pos  variant_list
#   mask_def.tsv   -- mask_name  annotation_categories
#
# Prerequisites:
#   1. Run create_ccre_bed.py locally to build blood_ccre_gene.bed
#   2. Upload blood_ccre_gene.bed to DNAnexus
#   3. Compile this workflow:  ./dx_compile.sh
#
# Submit:
#   dx run /workflows/ccre_anno/ccre_anno \
#       -istage-common.step2_chrom_manifest=file-XXX \
#       -istage-common.blood_ccre_gene_bed=file-YYY \
#       --destination /resources/ccre_anno/ -y --brief
#
# Use the three output file IDs in your gene-based RVAS config:
#   anno_file: file-AAA
#   set_list:  file-BBB
#   mask_def:  file-CCC

workflow ccre_anno {

    input {
        # Gene-based manifest: 7-column TSV (no header)
        # chrom[0], pvar_name[1], psam_name[2], pgen_name[3],
        # pvar_id[4], psam_id[5], pgen_id[6]
        File   step2_chrom_manifest

        # BED produced by create_ccre_bed.py (5 cols, no header):
        # chrom, start, end, gene, category
        File   blood_ccre_gene_bed

        # Variant ID prefix matching UKB pvar files (default: DRAGEN)
        String variant_id_prefix = "DRAGEN"
    }

    Array[Array[String]] chroms = read_tsv(step2_chrom_manifest)
    # Subtract 1 for the trailing empty row produced by read_tsv on a newline-terminated file
    Int n_chroms = length(chroms) - 1

    scatter(i in range(n_chroms)) {
        call annotate_chrom {
            input:
                pvar              = "dx://" + chroms[i][4],
                chrom_with_prefix = chroms[i][0],
                ccre_bed          = blood_ccre_gene_bed,
                variant_id_prefix = variant_id_prefix,
        }
    }

    call concatenate_anno {
        input:
            chrom_anno_files = annotate_chrom.chrom_anno_tsv,
    }

    output {
        File anno_file = concatenate_anno.anno_file
        File set_list  = concatenate_anno.set_list
        File mask_def  = concatenate_anno.mask_def
    }

}


# ---------------------------------------------------------------------------
# Task: annotate_chrom
#
# For one chromosome: stream the pvar file, filter to positions that fall
# within a cCRE interval, and write variant_id/gene/category rows.
# ---------------------------------------------------------------------------
task annotate_chrom {

    input {
        File   pvar
        File   ccre_bed
        String chrom_with_prefix
        String variant_id_prefix
    }

    command <<<
        echo "$(date -u) annotate_chrom ~{chrom_with_prefix} start"

        python3 << 'PYEOF'
import duckdb, sys

pvar_path         = "~{pvar}"
ccre_path         = "~{ccre_bed}"
chrom_with_prefix = "~{chrom_with_prefix}"
prefix            = "~{variant_id_prefix}"
out_path          = "chrom_" + chrom_with_prefix + "_ccre_anno.tsv"

con = duckdb.connect()

# pvar: lines starting with '#' are comments/header; positional columns
# column0=CHROM, column1=POS, column2=ID, column3=REF, column4=ALT
con.execute(f"""
    CREATE TABLE pvar AS
    SELECT
        '{prefix}' || ':' || column0 || ':' || CAST(column1 AS VARCHAR)
            || ':' || column3 || ':' || column4  AS variant_id,
        column0                                  AS chrom,
        CAST(column1 AS INTEGER)                 AS pos
    FROM read_csv(
        '{pvar_path}',
        delim   = '\t',
        comment = '#',
        header  = false,
        columns = {{
            'column0': 'VARCHAR',
            'column1': 'INTEGER',
            'column2': 'VARCHAR',
            'column3': 'VARCHAR',
            'column4': 'VARCHAR'
        }}
    )
""")

# Load only the cCRE intervals for this chromosome
con.execute(f"""
    CREATE TABLE ccre AS
    SELECT
        CAST(column1 AS INT) AS start,
        CAST(column2 AS INT) AS end,
        column3              AS gene,
        column4              AS category
    FROM read_csv(
        '{ccre_path}',
        delim   = '\t',
        header  = false,
        columns = {{
            'column0': 'VARCHAR',
            'column1': 'INTEGER',
            'column2': 'INTEGER',
            'column3': 'VARCHAR',
            'column4': 'VARCHAR'
        }}
    )
    WHERE column0 = '{chrom_with_prefix}'
""")

n_pvar = con.execute("SELECT COUNT(*) FROM pvar").fetchone()[0]
n_ccre = con.execute("SELECT COUNT(*) FROM ccre").fetchone()[0]
print(f"{chrom_with_prefix}: {n_pvar:,} pvar variants, {n_ccre:,} cCRE intervals",
      flush=True)

# Interval join: variant position within [start, end)
con.execute(f"""
    COPY (
        SELECT DISTINCT p.variant_id, c.gene, c.category
        FROM pvar p
        JOIN ccre c ON p.pos >= c.start AND p.pos < c.end
        ORDER BY c.gene, p.variant_id
    ) TO '{out_path}' (DELIMITER '\t', HEADER false)
""")

n_out = con.execute(
    f"SELECT COUNT(*) FROM read_csv('{out_path}', header=false, columns={{'column0':'VARCHAR','column1':'VARCHAR','column2':'VARCHAR'}})"
).fetchone()[0]
print(f"{chrom_with_prefix}: {n_out:,} variant-gene annotations -> {out_path}", flush=True)
PYEOF

        echo "$(date -u) annotate_chrom ~{chrom_with_prefix} done"
    >>>

    runtime {
        docker: "jweinstk/duckdb:master"
        dx_instance_type: "mem2_ssd1_v2_x2"
    }

    parameter_meta {
        pvar: {
            description: "input pvar",
            patterns: ["*.pvar"],
            stream: true
        }
        ccre_bed: {
            description: "cCRE BED file",
            patterns: ["*.bed"]
        }
    }

    output {
        File chrom_anno_tsv = "chrom_" + chrom_with_prefix + "_ccre_anno.tsv"
    }

}


# ---------------------------------------------------------------------------
# Task: concatenate_anno
#
# Merge per-chromosome annotation TSVs and write the three REGENIE gene-based
# annotation files: anno_file.tsv, set_list.tsv, mask_def.tsv.
# ---------------------------------------------------------------------------
task concatenate_anno {

    input {
        Array[File] chrom_anno_files
    }

    command <<<
        echo "$(date -u) concatenate_anno start"

        for f in ~{sep=' ' chrom_anno_files}; do
            ln -s "$f" .
        done

        python3 << 'PYEOF'
import duckdb

con = duckdb.connect()

# Union all per-chromosome annotation files
con.execute("""
    CREATE TABLE anno AS
    SELECT column0 AS variant_id, column1 AS gene, column2 AS category
    FROM read_csv(
        '*_ccre_anno.tsv',
        delim   = '\t',
        header  = false,
        columns = {'column0': 'VARCHAR', 'column1': 'VARCHAR', 'column2': 'VARCHAR'}
    )
""")

n = con.execute("SELECT COUNT(*) FROM anno").fetchone()[0]
n_genes = con.execute("SELECT COUNT(DISTINCT gene) FROM anno").fetchone()[0]
print(f"Loaded {n:,} variant-gene annotations across {n_genes:,} genes", flush=True)

# --- anno_file.tsv: variant_id  gene  category ---
con.execute("""
    COPY (
        SELECT variant_id, gene, category
        FROM anno
        ORDER BY gene, variant_id, category
    ) TO 'anno_file.tsv' (DELIMITER '\t', HEADER false)
""")
print(f"Written: anno_file.tsv  ({n:,} rows)", flush=True)

# --- set_list.tsv: gene  chrom  start_pos  variant_list ---
# Variant ID format: PREFIX:CHROM:POS:REF:ALT  (split_part is 1-based)
con.execute("""
    CREATE TABLE set_list_tbl AS
    SELECT
        gene,
        split_part(variant_id, ':', 2)                       AS chrom,
        MIN(CAST(split_part(variant_id, ':', 3) AS INTEGER)) AS start_pos,
        string_agg(DISTINCT variant_id, ',')                 AS variant_list
    FROM anno
    GROUP BY gene, split_part(variant_id, ':', 2)
""")
n_sets = con.execute("SELECT COUNT(*) FROM set_list_tbl").fetchone()[0]
con.execute("""
    COPY (
        SELECT gene, chrom, start_pos, variant_list
        FROM set_list_tbl
        ORDER BY chrom, start_pos
    ) TO 'set_list.tsv' (DELIMITER '\t', HEADER false)
""")
print(f"Written: set_list.tsv  ({n_sets:,} genes)", flush=True)

# --- mask_def.tsv: mask_name  annotation_categories ---
# Discover categories dynamically so this task works for any cCRE category set
cats = con.execute(
    "SELECT DISTINCT category FROM anno ORDER BY category"
).fetchall()
with open("mask_def.tsv", "w") as fh:
    for (cat,) in cats:
        fh.write(f"M_{cat}\t{cat}\n")
print(f"Written: mask_def.tsv  ({len(cats)} mask(s))", flush=True)
for (cat,) in cats:
    print(f"  M_{cat}  {cat}", flush=True)
PYEOF

        echo "$(date -u) concatenate_anno done"
    >>>

    runtime {
        docker: "jweinstk/duckdb:master"
        dx_instance_type: "mem2_ssd1_v2_x2"
    }

    output {
        File anno_file = "anno_file.tsv"
        File set_list  = "set_list.tsv"
        File mask_def  = "mask_def.tsv"
    }

}
