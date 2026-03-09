import sys
import yaml

config_file = sys.argv[1]

with open(config_file) as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)

cost_limit = cfg["cost_limit"]
project = cfg["project"]
priority = cfg["priority"]
step1_pvar = cfg["step1_pvar"]
step1_psam = cfg["step1_psam"]
step1_pgen = cfg["step1_pgen"]
step1_prefix = cfg["step1_prefix"]
covariates = cfg["covariates"]
phenotypes = cfg["phenotypes"]
covariate_string = cfg["covariate_string"]
categorical_covariate_string = cfg["categorical_covariate_string"]
final_folder = cfg["final_folder"]
concatenate = cfg["concatenate"]
fix_step2_header_for_rap = cfg["fix_step2_header_for_rap"]
workflow = cfg["workflow"]
step2_block_size = cfg["step2_block_size"]
minMAC = cfg["minMAC"]
step2_chunk_manifest = cfg["step2_chunk_manifest"]
plink2_binary = cfg["plink2_binary"]

if __name__ == "__main__":
    print(
        'uv run dx run {workflow} \
     -istage-common.step1_pvar={step1_pvar} \
     -istage-common.step1_pgen={step1_pgen} \
     -istage-common.step1_psam={step1_psam} \
     -istage-common.step1_prefix={step1_prefix} \
     -istage-common.covariate_string="{covariate_string}" \
     -istage-common.covariates="{covariates}" \
     -istage-common.categorical_covariate_string="{categorical_covariate_string}" \
     -istage-common.phenotypes="{phenotypes}" \
     -istage-common.concatenate_into_parquet="{concatenate}" \
     -istage-common.minMAC="{minMAC}" \
     -istage-common.step2_chunk_manifest="{step2_chunk_manifest}" \
     -istage-common.step2_block_size="{step2_block_size}" \
     -istage-common.fix_step2_header_for_rap="{fix_step2_header_for_rap}" \
     -istage-common.plink2_binary="{plink2_binary}" \
     --folder="{final_folder}" \
     --tag "regenie" \
     --priority {priority} \
     --cost-limit {cost} \
     -y \
     --brief'.format(
            workflow=workflow,
            step1_pvar=step1_pvar,
            step1_pgen=step1_pgen,
            step1_psam=step1_psam,
            step1_prefix=step1_prefix,
            covariate_string=covariate_string,
            categorical_covariate_string=categorical_covariate_string,
            covariates=covariates,
            phenotypes=phenotypes,
            minMAC=minMAC,
            step2_block_size=step2_block_size,
            step2_chunk_manifest=step2_chunk_manifest,
            concatenate=concatenate,
            fix_step2_header_for_rap=fix_step2_header_for_rap,
            plink2_binary=plink2_binary,
            final_folder=final_folder,
            priority=priority,
            cost=cost_limit,
        )
    )
