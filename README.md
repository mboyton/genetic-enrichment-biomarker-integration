# Genetic Enrichment and Biomarker Integration Pipeline

A modular, reproducible workflow for identifying genetically enriched cell populations and translating these findings into candidate therapeutic targets through transcriptional programme analysis and drug-target integration.

The pipeline was developed as part of a University of Bristol × RIKEN collaboration.

---

# Repository Structure

```
01_config.sh
01_genetic-enrichment-analysis.sh
01_genetic-enrichment-visualisation.R
02_biomarker_integration.R
README.md
```

---

# Workflow Overview

The workflow consists of two independent analytical stages.

## Stage 1 – Genetic Enrichment Analysis

**Scripts**

```
01_config.sh
01_genetic-enrichment-analysis.sh
01_genetic-enrichment-visualisation.R
```

Given:

- A CELLECT checkout (`CELLECT_DIR`)
- A pre-standardised GWAS summary statistics file
- A fixed global sample size (`FIXED_N`)
- A gene × cell-type specificity matrix (CELLEX ESµ)

the pipeline:

1. Validates required software, configuration and input files.
2. Ensures the LDSC munging environment is available.
3. Verifies the GWAS summary statistics format.
4. Munges GWAS summary statistics using `mtag_munge.py`.
5. Generates a run-specific CELLECT configuration file.
6. Executes the CELLECT-LDSC Snakemake workflow.
7. Produces a publication-ready visualisation of cell-type enrichment results.

---

## Stage 2 – Biomarker Integration

**Script**

```
02_biomarker_integration.R
```

Using the significant cell populations identified during Stage 1, the pipeline:

1. Loads transcriptional programme annotations.
2. Identifies significantly enriched immune cell populations.
3. Prioritises transcriptional hub genes using WGCNA connectivity (kME).
4. Queries Open Targets protein interaction data.
5. Identifies interacting proteins with approved therapeutics.
6. Filters high-confidence interactions.
7. Produces timestamped candidate drug target outputs together with run metadata for reproducibility.

---

# Stage 1 Input Requirements

## GWAS Summary Statistics

The GWAS summary statistics must:

- Be tab-delimited
- Contain the following header exactly:

```
SNP    A1    A2    BETA    SE    P    N    Z
```

- Be harmonisable to HapMap3 alleles
- Have an appropriate fixed global sample size (`FIXED_N`) specified within `01_config.sh`

This repository assumes all upstream GWAS preprocessing has already been completed.

---

## Cell-Type Specificity Matrix (CELLEX ESµ)

The specificity matrix must:

- Be comma-delimited
- Contain a header row
- Have gene IDs in the first column
- Use Ensembl gene identifiers (`ENSG...`)
- Exclude Ensembl version suffixes (e.g. `ENSG00000123456`, not `ENSG00000123456.12`)
- Contain no duplicate genes
- Show substantial overlap with the CELLECT reference gene coordinate file:

```
data/shared/gene_coordinates.GRCh37.ensembl_v91.txt
```

Expected structure:

```
Rows    = genes
Columns = cell populations / experimental conditions
Values  = cell-type specificity scores (ESµ)
```

---

# Stage 1 Configuration

The following parameters are intended to be edited between analyses:

```
CELLECT_DIR
SPEC_RAW
GWAS_IN

RUN_TAG
FIXED_N
SMK_CORES

SNAKEMAKE_BIN
CONDA_SH
MAMBA_BIN

CONDA_PREFIX_DIR
MUNGE_ENV
```

All analytical logic is contained within the pipeline scripts.

---

# Stage 1 Pre-flight Checklist

Before running CELLECT, verify the specificity matrix.

## 1. Confirm valid Ensembl Gene IDs

```bash
zcat specificity.csv.gz | tail -n +2 | head -n 2000 | cut -d',' -f1 | \
awk 'BEGIN{ens=0;tot=0} {tot++; if($0 ~ /^ENSG[0-9]+$/) ens++} END{print ens,"of",tot}'
```

All inspected rows should match.

---

## 2. Confirm no duplicated genes

```bash
zcat specificity.csv.gz | tail -n +2 | cut -d',' -f1 | sort | uniq -d | wc -l
```

Expected output:

```
0
```

---

## 3. Confirm overlap with CELLECT reference coordinates

```bash
COORD="data/shared/gene_coordinates.GRCh37.ensembl_v91.txt"

cut -f1 "$COORD" | tail -n +2 | sort -u > coord.txt
zcat specificity.csv.gz | tail -n +2 | cut -d',' -f1 | sort -u > spec.txt
comm -12 coord.txt spec.txt | wc -l
```

Substantial overlap (typically >80%) is expected.

---

# Genome Reference

The default coordinate reference used by the pipeline is:

```
GRCh37
Ensembl v91
```

If an alternative genome build is used, both the gene coordinate reference and specificity matrix must be generated using the same annotation release.

---

# Stage 2 Input Requirements

Stage 2 requires:

- CELLECT cell-type enrichment results generated during Stage 1
- A transcriptional programme database (RDS format)

The script automatically integrates these data with Open Targets protein interaction and approved drug annotations.

---

# Outputs

## Stage 1

Following completion, the pipeline generates:

- CELLECT prioritisation results
- Cell-type enrichment statistics
- Publication-ready visualisation
- Run-specific CELLECT configuration
- Complete execution log

Outputs are written under:

```
CELLECT-${RUN_TAG}/CELLECT-LDSC/
```

---

## Stage 2

Timestamped outputs are written to the configured output directory.

These include:

- High-confidence approved drug interaction table (`.csv`)
- Run metadata (`.rds`)

The metadata records:

- Input file locations
- Analysis parameters
- Statistical thresholds
- Numbers of significant cell populations
- Numbers of prioritised biomarkers
- Numbers of protein interactions
- Numbers of candidate drug targets