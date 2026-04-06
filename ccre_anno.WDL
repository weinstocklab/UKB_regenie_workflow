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
# Uses DuckDB CLI directly (no Python in this image).
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

        # Extract first 5 pvar columns (CHROM,POS,ID,REF,ALT) to a local CSV.
        # This bypasses DuckDB's sniffer and avoids reading the INFO column.
        awk '!/^#/{OFS=","; print $1,$2,$3,$4,$5}' ~{pvar} > pvar_5col.csv
        echo "pvar columns extracted"

        # Extract cCRE intervals for this chromosome, sorted by start position.
        # Result is tiny: typically <10k rows.
        awk -v chrom="~{chrom_with_prefix}" '$1==chrom{print $2"\t"$3"\t"$4"\t"$5}' \
            ~{ccre_bed} | sort -k1,1n > chr_ccre.bed
        n_ccre=$(wc -l < chr_ccre.bed)
        echo "n_ccre=${n_ccre}"

        # Sweep-line interval join: O(n_pvar + n_ccre) time, O(n_ccre) memory.
        # Loads chr_ccre.bed into arrays (NR==FNR pass), then streams pvar_5col.csv.
        # Both inputs are sorted by position; lo pointer skips expired intervals.
        # Outputs comma-separated variant_id,gene,category (no spaces in any field).
        awk -v prefix="~{variant_id_prefix}" -v chrom="~{chrom_with_prefix}" '
        BEGIN { lo = 1; n = 0 }
        NR==FNR {
            n++
            s[n] = $1 + 0
            e[n] = $2 + 0
            g[n] = $3
            c[n] = $4
            next
        }
        {
            split($0, f, ",")
            pos = f[2] + 0
            vid = prefix ":" chrom ":" f[2] ":" f[4] ":" f[5]
            for (i = lo; i <= n; i++) {
                if (s[i] > pos) break
                if (e[i] > pos) print vid "," g[i] "," c[i]
            }
            while (lo <= n && e[lo] <= pos) lo++
        }
        ' chr_ccre.bed pvar_5col.csv | sort -u > anno_rows.csv

        n_anno=$(wc -l < anno_rows.csv)
        echo "n_anno=${n_anno}"

        # Write anno_rows.csv to parquet via DuckDB (small file, no join needed).
        duckdb -c "
        COPY (
            SELECT column0 AS variant_id, column1 AS gene, column2 AS category
            FROM read_csv('anno_rows.csv', delim=',', header=false,
                          columns={'column0':'VARCHAR','column1':'VARCHAR','column2':'VARCHAR'})
        ) TO 'chrom_~{chrom_with_prefix}_ccre_anno.parquet' (FORMAT PARQUET)
        " || { echo "ERROR: parquet write failed for ~{chrom_with_prefix}"; exit 1; }

        # Verify the output parquet is valid.
        duckdb -c "SELECT COUNT(*) AS n_anno FROM parquet_scan('chrom_~{chrom_with_prefix}_ccre_anno.parquet')" \
            || { echo "ERROR: parquet validation failed for ~{chrom_with_prefix}"; exit 1; }

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
        File chrom_anno_tsv = "chrom_" + chrom_with_prefix + "_ccre_anno.parquet"
    }

}


# ---------------------------------------------------------------------------
# Task: concatenate_anno
#
# Merge per-chromosome annotation TSVs and write the three REGENIE gene-based
# annotation files: anno_file.tsv, set_list.tsv, mask_def.tsv.
# Uses DuckDB CLI directly (no Python in this image).
# ---------------------------------------------------------------------------
task concatenate_anno {

    input {
        Array[File] chrom_anno_files
    }

    command <<<
        echo "$(date -u) concatenate_anno start"

        # Build a DuckDB-format quoted file path list for read_parquet()
        FILES=""
        for f in ~{sep=' ' chrom_anno_files}; do
            [ -z "$FILES" ] && FILES="'$f'" || FILES="$FILES, '$f'"
        done

        mkdir -p /tmp/duckdb_spill

        duckdb << DUCKEOF
-- Cap DuckDB at 12 GB; spill to SSD rather than OOM.
SET memory_limit = '12GB';
SET temp_directory = '/tmp/duckdb_spill';

SELECT COUNT(*) AS n_rows, COUNT(DISTINCT gene) AS n_genes
FROM read_parquet([$FILES]);

-- anno_file.tsv: variant_id  gene  category
-- Read parquet directly (no intermediate table) to avoid double-materialising.
COPY (
    SELECT DISTINCT variant_id, gene, category
    FROM read_parquet([$FILES])
    ORDER BY gene, variant_id, category
) TO 'anno_file.tsv' (DELIMITER ' ', HEADER false);

-- set_list.tsv: gene  chrom  start_pos  variant_list
COPY (
    SELECT
        gene,
        split_part(variant_id, ':', 2)                       AS chrom,
        MIN(CAST(split_part(variant_id, ':', 3) AS INTEGER)) AS start_pos,
        string_agg(DISTINCT variant_id, ',')                 AS variant_list
    FROM read_parquet([$FILES])
    GROUP BY gene, split_part(variant_id, ':', 2)
    ORDER BY chrom, start_pos
) TO 'set_list.tsv' (DELIMITER ' ', HEADER false);

-- mask_def.tsv: mask_name  category (one mask per cell type, discovered dynamically)
COPY (
    SELECT 'M_' || category AS mask_name, category
    FROM (SELECT DISTINCT category FROM read_parquet([$FILES]) ORDER BY category)
) TO 'mask_def.tsv' (DELIMITER ' ', HEADER false);

SELECT mask_name, category FROM read_csv('mask_def.tsv', delim=' ', header=false,
    columns={'mask_name':'VARCHAR','category':'VARCHAR'});
DUCKEOF

        echo "$(date -u) concatenate_anno done"
    >>>

    runtime {
        docker: "jweinstk/duckdb:master"
        dx_instance_type: "mem2_ssd1_v2_x4"
    }

    output {
        File anno_file = "anno_file.tsv"
        File set_list  = "set_list.tsv"
        File mask_def  = "mask_def.tsv"
    }

}
