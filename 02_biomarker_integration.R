# ============================================================
# LOAD CELLECT, TRANSCRIPTIONAL PROGRAMS AND MR RESULTS
# ============================================================

suppressPackageStartupMessages({
 library(data.table)
 library(dplyr)
})

# ! paramteres #####
cellect_pval_threshold = 0.05
mr_method = "Inverse variance weighted"
run_output_directory = "C:/Users/uh21133/OneDrive - University of Bristol/Desktop/2026-06/Careers/Applications/NHS_Principal-Data-Scientist/interview_preparation_2/code_example/working"


# ============================================================
# FILE PATHS
# ============================================================

cellect_file_path <- paste0(
 "C:/Users/uh21133/OneDrive - University of Bristol/Desktop/",
 "2026-06/Careers/Applications/NHS_Principal-Data-Scientist/",
 "interview_preparation_2/code_example/working/",
 "LawlorPBMC__PGCMDD_Lawlor_20260219_1.cell_type_results.txt"
)

marker_database_path <- paste0(
 "C:/Users/uh21133/OneDrive - University of Bristol/Desktop/",
 "2026-06/Careers/Applications/NHS_Principal-Data-Scientist/",
 "interview_preparation_2/code_example/working/",
 "LawlorPBMC_wgcna_transcriptional_programs.rds"
)


# ============================================================
# CHECK FILES EXIST
# ============================================================

stopifnot(
 file.exists(cellect_file_path),
 file.exists(marker_database_path),
 file.exists(mr_file_path)
)

# ============================================================
# 1. LOAD CELLECT RESULTS
# ============================================================

cellect_results <- fread(
 cellect_file_path,
 data.table = FALSE
)

cellect_results <- cellect_results %>%
 mutate(
  Name = as.character(Name),
  Coefficient = as.numeric(Coefficient),
  Coefficient_std_error = as.numeric(Coefficient_std_error),
  Coefficient_P_value = as.numeric(Coefficient_P_value)
 )

cat(
 "\nCELLECT results loaded:",
 nrow(cellect_results),
 "cell types\n"
)

print(head(cellect_results))

# ============================================================
# 2. LOAD TRANSCRIPTIONAL PROGRAM DATABASE
# ============================================================

transcriptional_programs <- readRDS(
 marker_database_path
)

cat(
 "\nTranscriptional programmes loaded:",
 length(transcriptional_programs),
 "\n"
)

print(names(transcriptional_programs))


# ============================================================
# IDENTIFY SIGNIFICANT CELL TYPES AND EXTRACT TOP MARKERS
# ============================================================

top_n_markers <- 5

significant_cell_types <- cellect_results %>%
 filter(Coefficient_P_value < cellect_pval_threshold) %>%
 arrange(Coefficient_P_value)

cat(
 "\nSignificant CELLECT cell types:",
 nrow(significant_cell_types),
 "\n"
)

# Keep only significant cell types represented in the marker database
significant_programmes <- transcriptional_programs[
 intersect(
  significant_cell_types$Name,
  names(transcriptional_programs)
 )
]

# Combine programmes and select top WGCNA hub genes by kME
top_markers <- bind_rows(
 significant_programmes,
 .id = "cell_type"
) %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(!is.na(kME)) %>%
 group_by(cell_type) %>%
 arrange(desc(kME), .by_group = TRUE) %>%
 slice_head(n = top_n_markers) %>%
 ungroup()

print(
 top_markers %>%
  select(
   any_of(c(
    "cell_type",
    "gene",
    "module",
    "kME",
    "mean_logcpm_diff_stim_minus_base",
    "module_trait_fdr"
   ))
  )
)

# ============================================================
# BUILD PRIMARY TARGET + PROTEIN INTERACTOR TABLE
# ============================================================

library(dplyr)
library(purrr)
library(otargen)

primary_targets <- top_markers %>%
 distinct(gene, .keep_all = TRUE)

protein_interactions <- purrr::map_dfr(
 primary_targets$gene,
 function(primary_gene) {
  
  result <- interactionsQuery(
   ensgId = primary_gene,
   sourceDatabase = "intact",
   size = 25
  )
  
  if (is.null(result) || nrow(result) == 0) {
   return(NULL)
  }
  
  result %>%
   mutate(
    primary_target = primary_gene
   )
 }
)


names(protein_interactions)
head(protein_interactions)

# ============================================================
# QUERY KNOWN DRUGS FOR INTERACTING PROTEINS
# ============================================================

library(dplyr)
library(purrr)
library(tidyr)
library(otargen)

network_targets <- protein_interactions %>%
 mutate(
  network_gene = if_else(
   targetA.id == primary_target,
   targetB.id,
   targetA.id
  ),
  network_gene_symbol = if_else(
   targetA.id == primary_target,
   targetB.approvedSymbol,
   targetA.approvedSymbol
  )
 ) %>%
 filter(
  !is.na(network_gene),
  network_gene != primary_target
 ) %>%
 distinct(
  primary_target,
  network_gene,
  network_gene_symbol,
  .keep_all = TRUE
 )

network_drugs <- network_targets %>%
 distinct(network_gene, network_gene_symbol) %>%
 mutate(
  drug_results = map(
   network_gene,
   function(gene_id) {
    
    result <- knownDrugsGeneQuery(
     ensgId = gene_id
    )
    
    if (is.null(result) || nrow(result) == 0) {
     return(NULL)
    }
    
    result
   }
  )
 ) %>%
 filter(
  map_lgl(
   drug_results,
   ~ !is.null(.x)
  )
 ) %>%
 unnest(drug_results)

network_drugs

# ============================================================
# APPROVED DRUGS FOR NETWORK INTERACTORS
# ============================================================

approved_network_drugs <- network_targets %>%
 inner_join(
  network_drugs %>%
   filter(maxClinicalStage == "APPROVAL") %>%
   distinct(
    network_gene,
    network_gene_symbol,
    drug.id,
    drug.name,
    .keep_all = TRUE
   ),
  by = c("network_gene", "network_gene_symbol")
 ) %>%
 select(
  primary_target,
  network_gene,
  network_gene_symbol,
  drug.id,
  drug.name,
  maxClinicalStage,
  score
 ) %>%
 arrange(
  primary_target,
  network_gene_symbol,
  drug.name
 )

# ============================================================
# HIGH-CONFIDENCE APPROVED DRUGS (>0.8 interaction score)
# ============================================================

high_confidence_drugs <- approved_network_drugs %>%
 filter(score > 0.8) %>%
 distinct(
  primary_target,
  network_gene,
  network_gene_symbol,
  drug.id,
  drug.name,
  maxClinicalStage,
  score
 ) %>%
 arrange(
  desc(score),
  network_gene_symbol,
  drug.name
 )

high_confidence_drugs

# ============================================================
# SAVE OUTPUT
# ============================================================
write.csv(
 high_confidence_drug_names,
 file.path(
  run_output_directory,
  "high_confidence_approved_drug_names.csv"
 ),
 row.names = FALSE
)
