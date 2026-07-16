library(data.table)
library(ggplot2)

# ---- load data ----
args <- commandArgs(trailingOnly = TRUE)

cellect_results <- fread(args[1])
plot_title <- args[2]

# ---- helper: plot with 95% CI, solid/hollow by significance ----
plot_cellect_ci <- function(df, title, ylim_range = NULL,
                            p_cut = 0.05, adjust = FALSE, z = 1.96) {
 dt <- as.data.table(df)
 dt[, Name := gsub("^DICE__", "", Name)]
 
 # 95% CI
 dt[, `:=`(
  ci_lo = Coefficient - z * Coefficient_std_error,
  ci_hi = Coefficient + z * Coefficient_std_error
 )]
 
 # significance flag (nominal or FDR)
 if (adjust) {
  dt[, P_adj := p.adjust(Coefficient_P_value, method = "BH")]
  dt[, Sig := P_adj < p_cut]
  subtitle <- sprintf("%s (solid = FDR < %.02f)", title, p_cut)
 } else {
  dt[, Sig := Coefficient_P_value < p_cut]
  subtitle <- sprintf("%s (solid = P < %.02f)", title, p_cut)
 }
 
 # order by coefficient
 dt[, Name := factor(Name, levels = dt[order(Coefficient)]$Name)]
 
 p <- ggplot(dt, aes(x = Name, y = Coefficient)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0) +
  geom_point(
   data = dt[Sig == FALSE],
   shape = 21,
   size = 3,
   stroke = 1.1,
   fill = NA,
   color = "black"
  ) +
  geom_point(
   data = dt[Sig == TRUE],
   shape = 21,
   size = 3,
   stroke = 1.1,
   fill = "#4472C4",
   color = "black"
  ) +
  theme_minimal(base_size = 13) +
  labs(
   x = "Cell type",
   y = "CELLECT coefficient (heritability enrichment)",
   title = subtitle,
   caption = "Points = coefficient; whiskers = 95% CI; solid = significant"
  )
 
 if (is.null(ylim_range)) {
  p + coord_flip()
 } else {
  p + coord_flip(ylim = ylim_range)
 }
}

# ---- generate visualisation ----
cellect_plot <- plot_cellect_ci(
 cellect_results,
 plot_title
)

# ---- save visualisation ----
output_file <- file.path(
 dirname(args[1]),
 paste0(plot_title, "_genetic_enrichment.png")
)

ggsave(
 output_file,
 plot = cellect_plot,
 width = 10,
 height = 7,
 dpi = 300
)