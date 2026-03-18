
# Regenie WDL Workflow

This WDL workflow runs [regenie](https://rgcgithub.github.io/regenie/) for GWAS/RVAS. This has been tested on the UKB RAP, but not elsewhere.

This WDL parallelizes step 1 using the instructions [here](https://github.com/rgcgithub/regenie/wiki/Further-parallelization-for-level-0-models-in-Step-1). It calls 
regenie v4.1 via a docker container provided by the RGC developers.

## Available Workflows

### GWAS
- `regenie_qt.WDL`: Quantitative traits (linear regression)
- `regenie_bt.WDL`: Binary traits (logistic regression with Firth fallback)

Both GWAS workflows share the same structure and inputs, using a manifest file for step 2 chunks and plink2 filtering.

### Gene-based RVAS
- `regenie_gene_qt.WDL`: Quantitative traits
- `regenie_gene_bt.WDL`: Binary traits

Gene-based workflows run burden and variance-component tests (SKAT-O, ACATO) per chromosome using rare variant masks. Step 1 is identical to the GWAS workflows; step 2 filters to rare variants (MAF < 1%) per chromosome and runs REGENIE's `--anno-file`/`--set-list`/`--mask-def` pipeline. See **Gene-based annotation files** below for how to prepare masks.

### Annotation preprocessing
- `ccre_anno.WDL`: Intersects ENCODE rE2G cCRE enhancer-gene intervals with UKB WGS pvar files on DNAnexus to produce `anno_file`, `set_list`, and `mask_def` for cCRE-based blood RVAS. Run this once; the outputs plug directly into `regenie_gene_bt/qt.WDL`. See **cCRE annotation** below.

## Workflow Inputs

### GWAS (`regenie_qt.WDL`, `regenie_bt.WDL`)

Both workflows share the same input parameters:

- `File step1_pvar`: The pvar file for step 1. Presumably, based on array genotypes have been filtered for MAF, HWE, missingness and LD pruned. 
- `File step1_psam`: The psam file for step 1.
- `File step1_pgen`: The pgen file for step 1.
- `String step1_prefix`: The prefix for step 1 output files.
- `File step2_chunk_manifest`: The manifest file for step 2. This is a TSV listing chromosome regions and their corresponding pvar, psam, and pgen files. See below for details.
- `String covariate_string`: A string of covariate column names.
- `String categorical_covariate_string`: A string of categorical covariate column names.
- `File plink2_binary`: The plink2 binary file (zip archive). This is used to filter variants in step2 for improved performance.
- `File covariates`: The covariates file.
- `File phenotypes`: The phenotypes file.
- `Boolean concatenate_into_parquet`: Whether to concatenate the summary statistics into a single parquet file.
- `Int n_step1`: The number of "l0" jobs in step 1 to parallelize over.
- `Boolean fix_step2_header_for_rap`: The plink2 psam files from UKB RAP may have formatting issues. Setting this to true applies a workaround to fix the header.
- `Int minMAC`: The minimum minor allele count for variants in step 2 (default: 10).
- `Int threads`: Number of threads to use for computation (default: 8).
- `Int step1_block_size`: Block size for step 1 (default: 500).
- `Int step2_block_size`: Block size for step 2 (default: 250).


### Step2 chunk manifest

```
> glimpse(manifest)
Rows: 316
Columns: 10
$ chrom                  <chr> "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"…
$ start                  <dbl> 1, 10000001, 20000001, 30000001, 40000001, 50000001, 60000001, 70000001, 80000001, 90000001, 100000001, 110000001, 120000001, 13000000…
$ end                    <dbl> 10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 120000000, 130000000, …
$ range                  <glue> "1:1-10000000", "1:10000001-20000000", "1:20000001-30000000", "1:30000001-40000000", "1:40000001-50000000", "1:50000001-60000000", "1…
$ plink2_chrom_pvar_name <glue> "/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, PLINK format [500k release]/ukb24308_c1_b0_v1.pvar", "/Bulk/DRAGEN WGS/DRAGEN…
$ plink2_chrom_psam_name <glue> "/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, PLINK format [500k release]/ukb24308_c1_b0_v1.psam", "/Bulk/DRAGEN WGS/DRAGEN…
$ plink2_chrom_pgen_name <glue> "/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, PLINK format [500k release]/ukb24308_c1_b0_v1.pgen", "/Bulk/DRAGEN WGS/DRAGEN…
$ plink2_chrom_pvar      <chr> "file-Gjx9g7jJ8yxYff6BB47xf8f1", "file-Gjx9g7jJ8yxYff6BB47xf8f1", "file-Gjx9g7jJ8yxYff6BB47xf8f1", "file-Gjx9g7jJ8yxYff6BB47xf8f1", "f…
$ plink2_chrom_psam      <chr> "file-GzbyjvQJBx1VF1BGV57PGy64", "file-GzbyjvQJBx1VF1BGV57PGy64", "file-GzbyjvQJBx1VF1BGV57PGy64", "file-GzbyjvQJBx1VF1BGV57PGy64", "f…
$ plink2_chrom_pgen      <chr> "file-Gjx8kF0J8yxpF75qxKgg261X", "file-Gjx8kF0J8yxpF75qxKgg261X", "file-Gjx8kF0J8yxpF75qxKgg261X", "file-Gjx8kF0J8yxpF75qxKgg261X", "f…
```

### Step2 chrom manifest (gene-based)

The gene-based workflows use a simpler per-chromosome manifest (`step2_chrom_manifest`) with one row per chromosome and no genomic window coordinates:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome name (e.g. `chr1`) |
| `plink2_chrom_pvar_name` | pvar filename |
| `plink2_chrom_psam_name` | psam filename |
| `plink2_chrom_pgen_name` | pgen filename |
| `plink2_chrom_pvar` | DNAnexus file ID for pvar |
| `plink2_chrom_psam` | DNAnexus file ID for psam |
| `plink2_chrom_pgen` | DNAnexus file ID for pgen |

Typically 22–23 rows (autosomes + X). No header row. The same manifest is also used by `ccre_anno.WDL`.

### Gene-based RVAS (`regenie_gene_qt.WDL`, `regenie_gene_bt.WDL`)

Same step 1 inputs as GWAS, plus:

- `File step2_chrom_manifest`: TSV listing one chromosome per row (7 columns: `chrom`, `pvar_name`, `psam_name`, `pgen_name`, `pvar_id`, `psam_id`, `pgen_id`). One row per chromosome (22–23 rows), no start/end needed.
- `File anno_file`: Variant→gene+category annotation file (3 columns: `variant_id`, `gene`, `category`). Produced by `create_gene_masks.py` (AlphaMissense) or `ccre_anno.WDL` (cCRE).
- `File set_list`: Gene→variant list (4 columns: `gene`, `chrom`, `start_pos`, `variant_list`).
- `File mask_def`: Mask name→annotation categories (2 columns: `mask_name`, `categories`).
- `String aaf_bins`: AAF thresholds for burden masks (default: `"0.001,0.01"`).
- `String build_mask`: Mask aggregation — `"max"` (dominant-like) or `"sum"` (additive, default: `"max"`).
- `String vc_tests`: Variance component tests (default: `"skato,acato-full"`).
- `Float vc_max_aaf`: Max AAF for VC tests (default: `0.01`).
- `Float max_maf`: plink2 max MAF filter to retain rare variants (default: `0.01`).
- `Int min_mac`: plink2 min MAC filter (default: `1`).

#### Gene-based annotation files

**AlphaMissense masks** (missense variant functional scores):
```bash
uv run python create_gene_masks.py \
    --alphamissense AlphaMissense_hg38.tsv.gz \
    --uniprot-mapping uniprot_human_gene_names.tsv.gz \
    --thresholds 0.34 0.564 \
    --out-dir gene_masks/
```

**cCRE masks** (ENCODE enhancer-gene links, blood cell types):
```bash
# Download from ENCODE portal (no account required)
wget https://www.encodeproject.org/files/ENCFF950FTI/@@download/ENCFF950FTI.bed.gz  # K562
wget https://www.encodeproject.org/files/ENCFF055YXP/@@download/ENCFF055YXP.bed.gz  # GM12878
wget https://www.encodeproject.org/files/ENCFF300OCC/@@download/ENCFF300OCC.bed.gz  # CD34+

uv run python create_ccre_bed.py \
    --cell-type K562:ENCFF950FTI.bed.gz \
    --cell-type GM12878:ENCFF055YXP.bed.gz \
    --cell-type CD34:ENCFF300OCC.bed.gz \
    --score-threshold 0.5 \
    --out blood_ccre_gene.bed

dx upload blood_ccre_gene.bed --destination /resources/
# then run ccre_anno.WDL (see below)
```

### cCRE annotation (`ccre_anno.WDL`)

Intersects a cCRE BED file with UKB WGS pvar files on DNAnexus. Run once; outputs reuse across analyses.

- `File step2_chrom_manifest`: Same per-chromosome manifest used by the gene-based workflows.
- `File blood_ccre_gene_bed`: BED file produced by `create_ccre_bed.py` (5 columns: `chrom`, `start`, `end`, `gene`, `category`).
- `String variant_id_prefix`: Prefix for variant IDs matching pvar files (default: `"DRAGEN"`).

Outputs: `anno_file.tsv`, `set_list.tsv`, `mask_def.tsv` — pass their DNAnexus file IDs directly to `regenie_gene_bt/qt.WDL`.

```bash
./dx_compile.sh  # compile ccre_anno.WDL (one-time)
dx run /workflows/ccre_anno/ccre_anno \
    -istage-common.step2_chrom_manifest=file-XXX \
    -istage-common.blood_ccre_gene_bed=file-YYY \
    --destination /resources/ccre_anno/ -y --brief
```

## Workflow Outputs

- `File loco_list`: The list of LOCO predictions from step 1.
- `Array[File] locos`: The array of LOCO files from step 1.
- `Array[Array[File]] summary_stats`: The summary statistics files from step 2.
- `File parquet`: The concatenated files from step 2 into a single parquet file (optional).

## Dependencies

This project uses [uv](https://github.com/astral-sh/uv) for python dependencies (e.g., `dxpy`),
though this is not essential.

## Contact
Contact Josh Weinstock for details on the WDL. 
Note, this repo has no affiliation with the REGENIE developers; please contact them for questions about REGENIE itself.

If you use this, please credit the regenie developers. 
