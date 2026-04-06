version 1.0
workflow regenie_gene_bt {

    input {

        # Step 1 inputs (array genotypes for null model, same as GWAS workflow)
        File step1_pvar
        File step1_psam
        File step1_pgen
        String step1_prefix

        # Step 2 manifest (per-chromosome WGS pgen files)
        # TSV columns (0-based): chrom[0], pvar_name[1], psam_name[2], pgen_name[3],
        #                        pvar_id[4], psam_id[5], pgen_id[6]
        File step2_chrom_manifest

        String covariate_string
        String categorical_covariate_string
        File covariates
        File phenotypes

        # Gene-based test annotation files
        File anno_file    # variant -> gene + annotation category
        File set_list     # gene -> variant list
        File mask_def     # mask name -> annotation categories

        File plink2_binary

        # Step 1 parameters
        Int n_step1 = 16
        Int threads = 8
        Int step1_block_size = 500

        # Step 2 parameters
        Int step2_block_size = 250
        Int min_mac = 2          # min MAC for plink2 pre-filter (excludes singletons)
        Float max_maf = 0.01     # max MAF for plink2 pre-filter (retains rare variants)

        # Gene-based test parameters
        String aaf_bins = "0.001,0.01"       # AAF bins for burden masks
        String build_mask = "max"             # mask aggregation: max or sum
        String vc_tests = "skato,acato-full"  # variance component tests
        Float vc_max_aaf = 0.01              # max AAF for VC tests
        Boolean rgc_gene_p = true            # compute omnibus gene-level p-values

        Boolean concatenate_into_parquet = false
        Boolean fix_step2_header_for_rap = true
    }

    Array[Array[String]] step2_chroms = read_tsv(step2_chrom_manifest)
    Int n_chroms = length(step2_chroms) - 1

    call step1_split_l0 {

        input:
            pvar = step1_pvar,
            psam = step1_psam,
            pgen = step1_pgen,
            plink2_prefix = step1_prefix,
            covariates = covariates,
            covariate_string = covariate_string,
            categorical_covariate_string = categorical_covariate_string,
            phenotypes = phenotypes,
            bsize = step1_block_size,
            njobs = n_step1,
            threads = 4
    }

    scatter(j in range(n_step1)) {
        call step1_run_l0 {

            input:
                pvar = step1_pvar,
                psam = step1_psam,
                pgen = step1_pgen,
                plink2_prefix = step1_prefix,
                covariates = covariates,
                covariate_string = covariate_string,
                categorical_covariate_string = categorical_covariate_string,
                phenotypes = phenotypes,
                job = j + 1,
                master = step1_split_l0.master,
                snplist = step1_split_l0.snplist,
                bsize = step1_block_size,
                threads = threads
        }
    }

    call step1_run_l1 {
        input:
                pvar = step1_pvar,
                psam = step1_psam,
                pgen = step1_pgen,
                plink2_prefix = step1_prefix,
                covariates = covariates,
                covariate_string = covariate_string,
                categorical_covariate_string = categorical_covariate_string,
                phenotypes = phenotypes,
                master = step1_split_l0.master,
                snplist = step1_split_l0.snplist,
                level_0_list = flatten(step1_run_l0.level_0_list),
                bsize = step1_block_size,
                threads = threads
    }

    scatter(i in range(n_chroms)) {

        call step2 {

            input:
                pvar = "dx://" + step2_chroms[i][4],
                psam = "dx://" + step2_chroms[i][5],
                pgen = "dx://" + step2_chroms[i][6],
                plink2_prefix = basename(step2_chroms[i][1], ".pvar"),
                plink2_binary = plink2_binary,
                chrom = sub(step2_chroms[i][0], "chr", ""),
                out_prefix = "chrom_~{i + 1}",
                locos = step1_run_l1.locos,
                loco_list = step1_run_l1.loco_list,
                covariates = covariates,
                covariate_string = covariate_string,
                categorical_covariate_string = categorical_covariate_string,
                bsize = step2_block_size,
                min_mac = min_mac,
                max_maf = max_maf,
                threads = threads,
                fix_step2_header_for_rap = fix_step2_header_for_rap,
                phenotypes = phenotypes,
                anno_file = anno_file,
                set_list = set_list,
                mask_def = mask_def,
                aaf_bins = aaf_bins,
                build_mask = build_mask,
                vc_tests = vc_tests,
                vc_max_aaf = vc_max_aaf,
                rgc_gene_p = rgc_gene_p
        }
    }

    if(concatenate_into_parquet) {
        call concatenate {
            input:
                summary_stats = flatten(step2.rg_output)
        }
    }


    output {

        File loco_list                   = step1_run_l1.loco_list
        Array[File] locos                = step1_run_l1.locos
        Array[Array[File]] summary_stats = step2.rg_output
        Array[File] step2_log            = step2.log
        File? parquet                    = concatenate.parquet
        File step1_log                   = step1_run_l1.log

    }

}

task step1_split_l0 {

    input {

        File psam
        File pvar
        File pgen
        String plink2_prefix
        File covariates
        String covariate_string
        String categorical_covariate_string
        File phenotypes
        Int bsize
        Int threads
        Int njobs
    }

    command <<<

        echo "$(date -u) start now"

        ln -s ~{psam} ~{plink2_prefix}.psam
        ln -s ~{pvar} ~{plink2_prefix}.pvar
        ln -s ~{pgen} ~{plink2_prefix}.pgen

        regenie --step 1 \
          --bt \
          --loocv \
          --covarFile ~{covariates} \
          --covarColList ~{covariate_string} \
          --catCovarList ~{categorical_covariate_string} \
          --phenoFile ~{phenotypes} \
          --bsize ~{bsize} \
          --gz \
          --pgen ~{plink2_prefix} \
          --threads ~{threads} \
          --out fit_bin_l0 \
          --split-l0 fit_bin_parallel,~{njobs}

        ls -alh fit_bin_l0
        echo "$(date -u) done now"
    >>>

    runtime {
        docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        dx_instance_type: "mem2_ssd1_v2_x4"
    }

    output {
        File master = "fit_bin_parallel.master"
        Array[File] snplist = glob("fit_bin_parallel_job*.snplist")
    }

}

task step1_run_l0 {

    input {

        File psam
        File pvar
        File pgen
        String plink2_prefix
        File covariates
        String covariate_string
        String categorical_covariate_string
        File phenotypes
        File master
        Array[File] snplist
        Int job
        Int bsize
        Int threads
    }

    command <<<

        echo "$(date -u) start now"

        ln -s ~{psam} ~{plink2_prefix}.psam
        ln -s ~{pvar} ~{plink2_prefix}.pvar
        ln -s ~{pgen} ~{plink2_prefix}.pgen
        ln -s ~{master} fit_bin_parallel.master

        for list in ~{sep=' ' snplist}; do
            ln -s ${list} .
        done

        ls -alh

        regenie --step 1 \
          --bt \
          --loocv \
          --covarFile ~{covariates} \
          --covarColList ~{covariate_string} \
          --catCovarList ~{categorical_covariate_string} \
          --phenoFile ~{phenotypes} \
          --bsize ~{bsize} \
          --gz \
          --pgen ~{plink2_prefix} \
          --threads ~{threads} \
          --run-l0 fit_bin_parallel.master,~{job} \
          --out fit_bin_l0_~{job}

        ls -alh fit_bin_l0_~{job}
        echo "$(date -u) done now"
    >>>

    runtime {
        docker: "ghcr.io/rgcgithub/regenie/regenie:v4.0.gz"
        dx_instance_type: "mem2_ssd1_v2_x8"
    }

    output {
        Array[File] level_0_list = glob("fit_bin_parallel_job~{job}_l0_Y*")
    }

}

task step1_run_l1 {

    input {

        File psam
        File pvar
        File pgen
        String plink2_prefix
        File covariates
        String covariate_string
        String categorical_covariate_string
        File phenotypes
        File master
        Array[File] snplist
        Array[File] level_0_list
        Int bsize
        Int threads
    }

    command <<<

        echo "$(date -u) start now"

        ln -s ~{psam} ~{plink2_prefix}.psam
        ln -s ~{pvar} ~{plink2_prefix}.pvar
        ln -s ~{pgen} ~{plink2_prefix}.pgen
        ln -s ~{master} fit_bin_parallel.master

        for list in ~{sep=' ' snplist}; do
            ln -s ${list} .
        done

        for list in ~{sep=' ' level_0_list}; do
            echo $list
            ln -s ${list} .
        done

        regenie --step 1 \
          --bt \
          --loocv \
          --covarFile ~{covariates} \
          --covarColList ~{covariate_string} \
          --catCovarList ~{categorical_covariate_string} \
          --phenoFile ~{phenotypes} \
          --bsize ~{bsize} \
          --gz \
          --pgen ~{plink2_prefix} \
          --threads ~{threads} \
          --run-l1 fit_bin_parallel.master \
          --out fit_bin_l1

        # remove /home/dnanexus/work from path
        sed -i 's|/home/dnanexus/work/||g' fit_bin_l1_pred.list

        mv ~/work/*loco.gz .

        ls -alh

        echo "$(date -u) done now"
    >>>

    runtime {
        docker: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
        dx_instance_type: "mem2_ssd1_v2_x8"
    }

    output {
        Array[File] locos = glob("*loco.gz")
        File loco_list    = "fit_bin_l1_pred.list"
        File log          = "fit_bin_l1.log"
    }

}

task step2 {

    input {
        File psam
        File pvar
        File pgen
        String plink2_prefix
        String covariate_string
        String categorical_covariate_string
        File covariates
        File phenotypes
        File plink2_binary
        File anno_file
        File set_list
        File mask_def
        Array[File] locos
        File loco_list
        String out_prefix
        String chrom
        Int bsize
        Int min_mac
        Float max_maf
        Int threads
        Boolean fix_step2_header_for_rap
        String aaf_bins
        String build_mask
        String vc_tests
        Float vc_max_aaf
        Boolean rgc_gene_p
    }

    command <<<

        echo "$(date -u) Starting now"

        ln -s ~{psam} ~{plink2_prefix}.psam
        ln -s ~{pvar} ~{plink2_prefix}.pvar
        ln -s ~{pgen} ~{plink2_prefix}.pgen

        if [[ "~{fix_step2_header_for_rap}" == "true" ]]; then
            echo "$(date -u) Fixing psam header now for RAP files; here's the original psam:"
            head  ~{plink2_prefix}.psam

            awk 'NR==1 {print "#FID IID SEX"; print $1, $1, $5; next} {print $1, $1, $5}' ~{plink2_prefix}.psam > tmp.psam

            mv ~{plink2_prefix}.psam ~{plink2_prefix}.psam.old
            mv tmp.psam ~{plink2_prefix}.psam

            echo "$(date -u) Done fixing psam header now: "
            head  ~{plink2_prefix}.psam
        else
            echo "$(date -u) Fixing psam header now; here's the original psam:"
            head  ~{plink2_prefix}.psam

            awk 'NR==1 {print "#FID IID SEX"; next} {print $1, $1, $2}' ~{plink2_prefix}.psam > tmp.psam

            mv ~{plink2_prefix}.psam ~{plink2_prefix}.psam.old
            mv tmp.psam ~{plink2_prefix}.psam

            echo "$(date -u) Done fixing psam header now:"
            head  ~{plink2_prefix}.psam
        fi

        for loco in ~{sep=' ' locos}; do
            echo $loco
            ln -s ${loco} .
        done

        echo "$(date -u) Now filtering with plink2"

        apt-get update
        apt-get install unzip

        unzip ~{plink2_binary}

        # Keep rare variants for gene-based tests (no snps-only: retain indels/PTVs)
        ./plink2 --pfile ~{plink2_prefix} \
            --chr ~{chrom} \
            --max-maf ~{max_maf} \
            --mac ~{min_mac} \
            --make-pgen \
            --out ~{plink2_prefix}.filtered

        ls -alh
        echo "$(date -u) Done filtering with plink2; now running regenie gene-based tests"

        regenie \
          --step 2 \
          --pgen ~{plink2_prefix}.filtered \
          --covarFile ~{covariates} \
          --covarColList ~{covariate_string} \
          --catCovarList ~{categorical_covariate_string} \
          --phenoFile ~{phenotypes} \
          --bsize ~{bsize} \
          --bt \
          --firth --approx \
          --threads ~{threads} \
          --pred ~{loco_list} \
          --anno-file ~{anno_file} \
          --set-list ~{set_list} \
          --mask-def ~{mask_def} \
          --aaf-bins ~{aaf_bins} \
          --build-mask ~{build_mask} \
          --vc-tests ~{vc_tests} \
          --vc-maxAAF ~{vc_max_aaf} \
          --af-cc \
          ~{if rgc_gene_p then "--rgc-gene-p" else ""} \
          --gz \
          --out ~{out_prefix}

        ls -alh
        echo "$(date -u) Done now"
    >>>

    runtime {
        docker: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
        dx_instance_type: "mem2_ssd1_v2_x8"
    }

    parameter_meta {
        pgen: {
            description: "input pgen",
            patterns: ["*.pgen"],
            stream: true
        }
        pvar: {
            description: "input pvar",
            patterns: ["*.pvar"],
            stream: true
        }
        psam: {
            description: "input psam",
            patterns: ["*.psam"],
            stream: true
        }
    }

    output {
        Array[File] rg_output = glob("*.regenie.gz")
        File log = "~{out_prefix}.log"
    }

}

task concatenate {

    input {
        Array[File] summary_stats
    }

    command <<<

        for l in ~{sep=' ' summary_stats}; do
            ln -s ${l} .
        done

        duckdb -c "COPY
         (SELECT
             concat('chr', CHROM) AS CHROM,
             GENPOS AS POS,
             ID,
             ALLELE0,
             ALLELE1,
             A1FREQ,
             N,
             TEST,
             BETA,
             SE,
             CHISQ,
             LOG10P,
             EXTRA,
             replace(regexp_replace(filename, '^.*chrom_[0-9]+_', ''), '.regenie.gz', '') AS phenotype
             FROM read_csv('*.gz',
                sep = ' ',
                comment = '#',
                filename = true,
                nullstr = 'NA',
                union_by_name = true,
                types={'ALLELE0': 'VARCHAR', 'ALLELE1': 'VARCHAR', 'ID': 'VARCHAR', 'EXTRA': 'VARCHAR',
                       'A1FREQ': 'DOUBLE', 'BETA': 'DOUBLE', 'SE': 'DOUBLE'},
                ignore_errors = true
            )
         )
         TO 'concatenated_results.parquet' (FORMAT PARQUET) "
    >>>

    runtime {
        docker: "jweinstk/duckdb:master"
        dx_instance_type: "mem2_ssd1_v2_x2"
    }

    output {
        File parquet = "concatenated_results.parquet"
    }

}
