# ============================================================
# LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({
 library(data.table)
 library(dplyr)
 library(purrr)
 library(tidyr)
 library(otargen)
})

# ============================================================
# CONFIGURATION
# ============================================================

cellect_pval_threshold = 0.05
top_n_markers = 5
drug_confidence_threshold = 0.8

run_output_directory = "/path/to/output_directory"

cellect_file_path = "/path/to/cellect_results.cell_type_results.txt"

marker_database_path = "/path/to/transcriptional_programs.rds"


# ============================================================
# CHECK FILES EXIST
# ============================================================

stopifnot(
 file.exists(cellect_file_path),
 file.exists(marker_database_path)
)

dir.create(
 run_output_directory,
 recursive = TRUE,
 showWarnings = FALSE
)

# ============================================================
# 1. Load genetic enrichment results
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

# ============================================================
# 2. Load biomarker database
# ============================================================

transcriptional_programs <- readRDS(
 marker_database_path
)

cat(
 "\nTranscriptional programmes loaded:",
 length(transcriptional_programs),
 "\n"
)

# ============================================================
# 3. Identify significant cell populations and markers
# ============================================================

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

# Combine programmes and select top hub genes by kME
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
# 4. Build primary target and protein interactor data
# ============================================================

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

cat(
 "\nProtein interactions identified:",
 nrow(protein_interactions),
 "\n"
)

# ============================================================
# 5. Query known drugs for interacting proteins (via API)
# ============================================================

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

cat(
 "\nDrug records identified:",
 nrow(network_drugs),
 "\n"
)

# ============================================================
# 6. Filter for approved drugs
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
# 7. Filter for high confidence, unique approved drugs
# ============================================================

high_confidence_drugs <- approved_network_drugs %>%
 filter(score > drug_confidence_threshold) %>%
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

cat(
 "\nHigh-confidence approved drug interactions:",
 nrow(high_confidence_drugs),
 "\n"
)


# ============================================================
# 8. Save outputs and run metadata
# ============================================================

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

output_file <- file.path(
 run_output_directory,
 paste0(
  "high_confidence_approved_drugs_",
  run_timestamp,
  ".csv"
 )
)

write.csv(
 high_confidence_drugs,
 output_file,
 row.names = FALSE
)

run_metadata <- list(
 run_timestamp = run_timestamp,
 cellect_file_path = cellect_file_path,
 marker_database_path = marker_database_path,
 cellect_pval_threshold = cellect_pval_threshold,
 top_n_markers = top_n_markers,
 drug_confidence_threshold = drug_confidence_threshold,
 significant_cell_types = nrow(significant_cell_types),
 primary_targets = nrow(primary_targets),
 protein_interactions = nrow(protein_interactions),
 high_confidence_drugs = nrow(high_confidence_drugs)
)

metadata_file <- file.path(
 run_output_directory,
 paste0(
  "biomarker_integration_metadata_",
  run_timestamp,
  ".rds"
 )
)

saveRDS(
 run_metadata,
 metadata_file
)

cat(
 "\nOutputs saved:",
 "\n -", output_file,
 "\n -", metadata_file,
 "\n"
)
