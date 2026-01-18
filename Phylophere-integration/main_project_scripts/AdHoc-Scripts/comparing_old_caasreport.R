
## 3. Overlap Analysis between Old and New CAAS Data

```{r caas_overlaps, fig.width=10, fig.height=10, dev="png", dpi=320}
# Old CAAS dataset (pre-data incident)
aggCAAS_df.old <- list.files(
  file.path(resultsDir, paste0("results.bak.16052025.pre-data-drama/2.CAAS/4.Randomizations/"), trait),
  pattern = "\\_caas.csv$", full.names = TRUE
) %>%
  lapply(read_csv) %>%
  bind_rows() %>%
  distinct() %>%
  mutate(GenePos = paste0(gene, "_", msa_pos))

# Background gene set
gene_bg <- read.table(file.path(dataDir, "2.Alignments/ordered_genes.txt"))$V1
bg_gene_n <- length(gene_bg)

# Define group subsets for both datasets
define_groups <- function(df) {
  list(
    ALL              = df,
    SC       = df %>% filter((isAntiCAAS == FALSE | isAntiCAAS == "CONS_ANTICAAS") & isSignificant),
    SC_PAT1  = df %>% filter((isAntiCAAS == FALSE | isAntiCAAS == "CONS_ANTICAAS") & isSignificant & pattern == "pattern1")
  )
}

groups.new <- define_groups(aggCAAS_df)
groups.old <- define_groups(aggCAAS_df.old)

# Hypergeometric test
hypergeom_test <- function(overlap, total_new, total_old, background) {
  phyper(q = overlap - 1, 
         m = total_old, 
         n = background - total_old, 
         k = total_new, 
         lower.tail = FALSE)
}

# Overlap analysis
analyze_overlap_detailed <- function(group_name) {
  new_df <- groups.new[[group_name]]
  old_df <- groups.old[[group_name]]
  
  new_genes <- unique(new_df$gene)
  old_genes <- unique(old_df$gene)
  overlap_genes <- intersect(new_genes, old_genes)
  only_new_genes <- setdiff(new_genes, old_genes)
  only_old_genes <- setdiff(old_genes, new_genes)
  
  new_caas <- unique(new_df$GenePos)
  old_caas <- unique(old_df$GenePos)
  overlap_caas <- intersect(new_caas, old_caas)
  only_new_caas <- setdiff(new_caas, old_caas)
  only_old_caas <- setdiff(old_caas, new_caas)
  
  tibble(
    Group = group_name,
    Category = c("ONLY_OLD", "ONLY_NEW", "OVERLAP"),
    Genes = c(length(only_old_genes), length(only_new_genes), length(overlap_genes)),
    CAAS  = c(length(only_old_caas),  length(only_new_caas),  length(overlap_caas))
  )
}

analyze_overlap_summary <- function(group_name) {
  new_df <- groups.new[[group_name]]
  old_df <- groups.old[[group_name]]
  
  # Unique CAAS entries
  new_caas <- new_df %>% distinct(GenePos, .keep_all = TRUE)
  old_caas <- old_df %>% distinct(GenePos, .keep_all = TRUE)
  
  # Overlapping CAAS (shared GenePos)
  overlap_gene_pos <- intersect(new_caas$GenePos, old_caas$GenePos)
  
  overlap_genes <- new_df %>%
    filter(GenePos %in% overlap_gene_pos) %>%
    pull(gene) %>%
    unique()
  
  n_overlap_genes <- length(overlap_genes)
  n_overlap_caas  <- length(overlap_gene_pos)
  
  # Only-NEW CAAS (non-overlapping)
  only_new_caas <- new_caas %>%
    filter(!(GenePos %in% overlap_gene_pos)) %>%
    rename(Gene = gene)
  
  # Only-OLD CAAS (non-overlapping)
  only_old_caas <- old_caas %>%
    filter(!(GenePos %in% overlap_gene_pos)) %>%
    rename(Gene = gene)
  
  # Exclude genes that appear in both only_new and only_old CAAS
  intersecting_genes <- intersect(only_new_caas$Gene, only_old_caas$Gene)
  intersecting_caas <- intersect(only_new_caas$GenePos, only_old_caas$GenePos)
  
  only_new_caas_filtered <- only_new_caas %>%
    filter(!(Gene %in% intersecting_genes))
  
  only_old_caas_filtered <- only_old_caas %>%
    filter(!(Gene %in% intersecting_genes))
  
  # Count genes and CAAS
  n_only_new_genes <- length(unique(only_new_caas_filtered$Gene))
  n_only_new_caas  <- nrow(only_new_caas)
  
  n_only_old_genes <- length(unique(only_old_caas_filtered$Gene))
  n_only_old_caas  <- nrow(only_old_caas)
  
  # Final tibble
  res <- tibble(
    Group = group_name,
    Category = c("ONLY_NEW", "ONLY_OLD", "OVERLAP"),
    Genes = c(n_only_new_genes, n_only_old_genes, length(intersecting_genes)),
    CAAS  = c(n_only_new_caas,  n_only_old_caas,  length(intersecting_caas))
  )
  
  return(res)
}

# Run overlap analysis for target groups
group_names <- c("ALL", "SC", "SC_PAT1")
detailed_overlap <- map_dfr(group_names, analyze_overlap_detailed)
overlap_summary <- map_dfr(group_names, analyze_overlap_summary)

# Display table
print(detailed_overlap)
print(overlap_summary)

# Save summary
outdir <- file.path(resultsDir, paste0("2.CAAS/5.Overlaps_old-new/", trait))
createDir(outdir)
write_csv(detailed_overlap, file.path(outdir, paste0(trait, "detailed_overlap.csv")))

# Plot stacked bars
plot_stacked_bar_percent <- function(df, var, title_suffix) {
  df_percent <- df %>%
    group_by(Group) %>%
    mutate(
      Total = sum(.data[[var]]),
      Percent = .data[[var]] / Total * 100,
      Label = .data[[var]]
    )
  
  ggplot(df_percent, aes(x = Group, y = Percent, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Label),
              position = position_stack(vjust = 0.5),
              color = "white", size = 3.5, fontface = "bold") +
    labs(
      title = paste("Proportional Overlap of", title_suffix, "-", trait),
      x = "Subset Group", y = paste("Proportion of", title_suffix, "(%)"), fill = "Category"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 100, by = 20)) +
    scale_fill_manual(values = c("ONLY_OLD" = "salmon3", "ONLY_NEW" = "skyblue3", "OVERLAP" = "forestgreen")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

# Plot and save
p_genes <- plot_stacked_bar_percent(detailed_overlap, "Genes", "Genes (OLD v NEW)")
p_caas  <- plot_stacked_bar_percent(detailed_overlap, "CAAS", "CAAS (OLD v NEW)")

p <- p_genes + p_caas

p

p_genes2 <- plot_stacked_bar_percent(overlap_summary, "Genes", "Genes (non-overlap)")
p_caas2  <- plot_stacked_bar_percent(overlap_summary, "CAAS", "CAAS (non-overlap)")

p2 <- p_genes2 + p_caas2

p2

ggsave(file.path(outdir, paste0(trait, "_genes_stacked_barplot.png")), p_genes, width = 7, height = 7, dpi = 320)
ggsave(file.path(outdir, paste0(trait, "_caas_stacked_barplot.png")), p_caas,  width = 7, height = 7, dpi = 320)
ggsave(file.path(outdir, paste0(trait, "_genes_caas_stacked_barplot.png")), p, width = 12, height = 7, dpi = 320)

ggsave(file.path(outdir, paste0(trait, "_genes_stacked_barplot2.png")), p_genes2, width = 7, height = 7, dpi = 320)
ggsave(file.path(outdir, paste0(trait, "_caas_stacked_barplot2.png")), p_caas2,  width = 7, height = 7, dpi = 320)
ggsave(file.path(outdir, paste0(trait, "_genes_caas_stacked_barplot2.png")), p2, width = 12, height = 7, dpi = 320)

```
