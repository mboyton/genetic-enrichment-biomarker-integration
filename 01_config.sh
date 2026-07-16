########################################
# RUN-SPECIFIC SETTINGS
########################################

CELLECT_DIR="/path/to/CELLECT"
SPEC_RAW="/path/to/cell_specificity_file.esmu.csv.gz"
GWAS_IN="/path/to/gwas_input.tsv"

RUN_TAG="example_run"
FIXED_N=50000 # Effective/sample size used for LDSC; amend per GWAS
SMK_CORES=4

########################################
# SOFTWARE / ENVIRONMENT PATHS
########################################

SNAKEMAKE_BIN="/path/to/snakemake"
CONDA_SH="/path/to/conda.sh"
MAMBA_BIN="/path/to/mamba"

CONDA_PREFIX_DIR="/path/to/snakemake_envs"
MUNGE_ENV="/path/to/munge_ldsc_env"