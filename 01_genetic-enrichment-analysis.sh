```bash
#!/bin/bash

########################################
# LOAD CONFIGURATION
########################################

set -euo pipefail

CONFIG_FILE="${1:-01_config.sh}"

[ -f "$CONFIG_FILE" ] || {
  echo "ERROR: Configuration file not found: $CONFIG_FILE" >&2
  exit 1
}

source "$CONFIG_FILE"

########################################
# DEFINE RUN-SPECIFIC PATHS
########################################

BASE_OUT="$CELLECT_DIR/CELLECT-${RUN_TAG}"
LDSC_OUT="$BASE_OUT/CELLECT-LDSC"

# Store run-specific munged GWAS inputs under CELLECT/custom_inputs.
# The cell-specificity file is referenced in place rather than copied.
CUSTOM_IN="$CELLECT_DIR/custom_inputs/${RUN_TAG}"

mkdir -p "$CUSTOM_IN" "$BASE_OUT"

########################################
# DEFINE OUTPUT FILES
########################################

# LDSC appends .sumstats.gz to this prefix.
MUNGED_PREFIX="$CUSTOM_IN/${RUN_TAG}"

CONFIG_YAML="$CELLECT_DIR/config_${RUN_TAG}.yml"

RUN_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$CELLECT_DIR/cellect_ldsc_run_${RUN_TAG}_${RUN_TIMESTAMP}.log"

########################################
# PRE-RUN CHECKS
########################################
[ -d "$CELLECT_DIR" ] || { echo "ERR: CELLECT_DIR not found: $CELLECT_DIR"; exit 1; }
[ -f "$SPEC_RAW" ] || { echo "ERR: ESmu file not found: $SPEC_RAW"; exit 1; }
[ -x "$SNAKEMAKE_BIN" ] || { echo "ERR: snakemake not found at $SNAKEMAKE_BIN"; exit 1; }
[ -f "$CONDA_SH" ] || { echo "ERR: conda shim not found at $CONDA_SH"; exit 1; }
[ -x "$MAMBA_BIN" ] || { echo "ERR: mamba not found at $MAMBA_BIN"; exit 1; }
[ -f "$GWAS_IN" ] || { echo "ERR: GWAS input not found: $GWAS_IN"; exit 1; }

########################################
# 1) Load conda
########################################
source "$CONDA_SH"

########################################
# 2) Ensure the LDSC munge env exists
########################################
if [ ! -d "$MUNGE_ENV" ]; then
  echo "[INFO] Creating LDSC munge env at: $MUNGE_ENV"
  conda env create -f "$CELLECT_DIR/ldsc/environment_munge_ldsc.yml" -p "$MUNGE_ENV"
fi
echo "[OK] munge env ready: $MUNGE_ENV"

########################################
# 3) Check GWAS header
########################################
echo "[INFO] Checking pre-standardised GWAS header in: $GWAS_IN"

awk 'BEGIN{FS=OFS="\t"}
NR==1{
  expected = "SNP\tA1\tA2\tBETA\tSE\tP\tN\tZ"
  got = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8
  if (got != expected) {
    print "ERR: Header mismatch." > "/dev/stderr"
    print "ERR: Expected: " expected > "/dev/stderr"
    print "ERR: Got:      " got      > "/dev/stderr"
    exit 1
  }
  exit 0
}' "$GWAS_IN"

echo "[OK] GWAS header looks good; using file in-place."

# Use the GWAS file directly for LDSC munging
GWAS_FOR_MUNGE="$GWAS_IN"

########################################
# 4) GWAS already has Z; use it for signed-sumstats
########################################
SIGNED_OPT=( --signed-sumstats Z,0 )

########################################
# 5) Always use a fixed global N for LDSC
########################################
if [ -z "$FIXED_N" ]; then
  echo "ERR: FIXED_N is empty but a global N is required."
  exit 1
fi
echo "[OK] Using fixed global N = $FIXED_N for LDSC"

########################################
# 6) Munge with LDSC (mtag_munge.py)
########################################
echo "[INFO] Munging GWAS -> ${MUNGED_PREFIX}.sumstats.gz"

conda run -p "$MUNGE_ENV" python "$CELLECT_DIR/ldsc/mtag_munge.py" \
  --sumstats "$GWAS_FOR_MUNGE" \
  --snp SNP --a1 A1 --a2 A2 --p P \
  --merge-alleles "$CELLECT_DIR/data/ldsc/w_hm3.snplist" \
  --keep-pval \
  --out "$MUNGED_PREFIX" \
  --n-value "$FIXED_N" \
  --signed-sumstats Z,0

[ -s "${MUNGED_PREFIX}.sumstats.gz" ] || { echo "ERR: munged sumstats not produced"; exit 1; }
echo "[OK] Munged sumstats: ${MUNGED_PREFIX}.sumstats.gz"

########################################
# 7) Write config.yml for this run
########################################
echo "[INFO] Writing config -> $CONFIG_YAML"
cat > "$CONFIG_YAML" <<YML
---
BASE_OUTPUT_DIR: ${BASE_OUT}

SPECIFICITY_INPUT:
  - id: cellect_output
    path: ${SPEC_RAW}

GWAS_SUMSTATS:
  - id: ${RUN_TAG}
    path: ${MUNGED_PREFIX}.sumstats.gz

ANALYSIS_TYPE:
  prioritization: True
  conditional: False
  heritability: False
  heritability_intervals: False

WINDOW_DEFINITION:
  WINDOW_SIZE_KB:
    100

GENE_COORD_FILE:
  data/shared/gene_coordinates.GRCh37.ensembl_v91.txt

KEEP_ANNOTS:
  False

LDSC_CONST:
  DATA_DIR:
    data/ldsc
  LDSC_DIR:
    ldsc
  NUMPY_CORES:
    1
YML

########################################
# 8) Make sure directory is unlocked & clean start
########################################
cd "$CELLECT_DIR"
"$SNAKEMAKE_BIN" --unlock -s cellect-ldsc.snakefile --configfile "$CONFIG_YAML" || true
rm -rf ".snakemake/locks"

########################################
# 9) Run Snakemake (foreground — correct for SLURM)
########################################

echo "[INFO] Starting Snakemake. Log: $LOG_FILE"

"$SNAKEMAKE_BIN" \
  --use-conda \
  --conda-prefix "$CONDA_PREFIX_DIR" \
  --conda-frontend mamba \
  -j "$SMK_CORES" \
  -s cellect-ldsc.snakefile \
  --configfile "$CONFIG_YAML" \
  --rerun-incomplete \
  > "$LOG_FILE" 2>&1

########################################
# 10) Plot CELLECT results
########################################

CELLECT_RESULTS="$LDSC_OUT/out/prioritization/cellect_output__${RUN_TAG}.cell_type_results.txt"
CELLECT_PLOT="$LDSC_OUT/out/prioritization/${RUN_TAG}_genetic_enrichment.png"

PLOT_SCRIPT="$(dirname "$0")/01_genetic-enrichment-visualisation.R"

[ -f "$CELLECT_RESULTS" ] || {
  echo "ERR: CELLECT results not found: $CELLECT_RESULTS"
  exit 1
}

[ -f "$PLOT_SCRIPT" ] || {
  echo "ERR: Plotting script not found: $PLOT_SCRIPT"
  exit 1
}

echo "[INFO] Plotting CELLECT results"

Rscript "$PLOT_SCRIPT" \
  "$CELLECT_RESULTS" \
  "$RUN_TAG"

########################################
# PIPELINE COMPLETE
########################################

echo
echo "========================================"
echo " Genetic enrichment analysis complete"
echo "========================================"
echo
echo "Log file:"
echo "  $LOG_FILE"
echo
echo "Key outputs:"
echo "  - Prioritization summary: $LDSC_OUT/results/prioritization.csv"
echo "  - Cell type results:      $CELLECT_RESULTS"
echo "  - Visualisation:          $CELLECT_PLOT"
echo