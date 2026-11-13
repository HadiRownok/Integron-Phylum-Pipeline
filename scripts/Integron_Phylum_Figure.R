# ==============================================================================
# INTEGRON DISTRIBUTION BY NCBI PHYLUM (Figure 2D style)
# Input : Integron_Phylum_Association.xlsx (from Integron_Phylum_Pipeline.sh)
# Output: PDF (vector) + 600 DPI PNG
# Taxonomy = Kraken2 classification of each integron-carrying contig
# ==============================================================================

# Run in RStudio (edit file_path below) or from the command line:
#   Rscript Integron_Phylum_Figure.R <Integron_Phylum_Association.xlsx|.tsv> [out_dir] [n_samples]

# 1. Paths
args <- commandArgs(trailingOnly = TRUE)
file_path  <- if (length(args) >= 1) args[1] else
  "Integron_Phylum_results/Integron_Phylum_Association.xlsx"
output_dir <- if (length(args) >= 2) args[2] else dirname(file_path)
output_dir <- paste0(sub("/$", "", output_dir), "/")
is_xlsx    <- grepl("\\.xlsx$", file_path, ignore.case = TRUE)

# 2. Packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", if (is_xlsx) "readxl")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  if (!interactive()) stop("Missing R packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  if (is_xlsx) library(readxl)
})

# 3. Settings
min_percent       <- 1      # phyla below this % of integrons -> "Other"
show_unclassified <- TRUE   # FALSE = drop Unclassified integrons from the figure

# 4. Read
df <- if (is_xlsx) read_excel(file_path, sheet = "Integron_Phylum_Association") else
  read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!show_unclassified) df <- filter(df, Phylum != "Unclassified")

# 5. Group small phyla into "Other (n phyla)"
ph_pct <- df %>% count(Phylum) %>% mutate(pct = 100 * n / sum(n))
minor  <- ph_pct$Phylum[ph_pct$pct < min_percent & ph_pct$Phylum != "Unclassified"]
if (length(minor) > 1) {
  df$Phylum[df$Phylum %in% minor] <- paste0("Other (", length(minor), " phyla)")
}

# 6. Counts per phylum x integron type
summary_df <- df %>%
  count(Phylum, integron_type, name = "count") %>%
  complete(Phylum, integron_type = c("CALIN", "In0", "complete"), fill = list(count = 0))

phylum_totals <- summary_df %>%
  group_by(Phylum) %>%
  summarise(total = sum(count),
            label_text = paste(count[integron_type == "complete"],
                               count[integron_type == "In0"],
                               count[integron_type == "CALIN"], sep = " / "),
            .groups = "drop") %>%
  mutate(pct = round(100 * total / sum(total), 1),
         right_box_text = paste0(total, " (", pct, "%)"))

# samples in which each phylum carries >= 1 integron (cf. "genomes containing integrons" in Neron et al. 2022)
n_samples_all <- if (length(args) >= 3) as.integer(args[3]) else if (is_xlsx)
  tryCatch(nrow(read_excel(file_path, sheet = "Per_sample")), error = function(e) n_distinct(df$sample)) else
  n_distinct(df$sample)
samp <- df %>% group_by(Phylum) %>% summarise(n_samp = n_distinct(sample), .groups = "drop")
phylum_totals <- phylum_totals %>% left_join(samp, by = "Phylum") %>%
  mutate(samp_text = sprintf("%d/%d (%.0f%%)", n_samp, n_samples_all, 100 * n_samp / n_samples_all))

# y order (bottom -> top): Unclassified, Other, then phyla by size (largest on top)
special  <- phylum_totals$Phylum[phylum_totals$Phylum == "Unclassified" | grepl("^Other", phylum_totals$Phylum)]
special  <- special[order(special != "Unclassified")]
regular  <- phylum_totals %>% filter(!Phylum %in% special) %>% arrange(total) %>% pull(Phylum)
lv <- c(special, regular)

summary_df$Phylum    <- factor(summary_df$Phylum, levels = lv)
phylum_totals$Phylum <- factor(phylum_totals$Phylum, levels = lv)
summary_df$integron_type <- factor(summary_df$integron_type, levels = c("CALIN", "In0", "complete"))

# 7. Colours
type_colors <- c("complete" = "#3182bd", "In0" = "#de2d26", "CALIN" = "#fec44f")
xmax <- max(phylum_totals$total) * 1.35

# 8. Plot 1: stacked bars
p1 <- ggplot(summary_df, aes(y = Phylum, x = count, fill = integron_type)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.55, color = "black", linewidth = 0.3) +
  geom_text(data = phylum_totals, aes(y = Phylum, x = total + xmax * 0.02, label = label_text),
            inherit.aes = FALSE, size = 3.6, fontface = "bold", hjust = 0, color = "#1a1a1a") +
  scale_fill_manual(values = type_colors, name = "Integron Type", breaks = c("complete", "In0", "CALIN")) +
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  labs(title = "Integron Distribution by NCBI Phylum", x = "Number of Integrons", y = NULL) +
  theme_classic() +
  theme(
    plot.title   = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.text.y  = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
    legend.position   = "right",
    legend.background = element_rect(fill = "#f8f9fa", color = "#cccccc", linewidth = 0.4),
    legend.title = element_text(face = "bold", size = 10),
    legend.text  = element_text(size = 9.5),
    plot.margin  = margin(5, 15, 5, 5))

# 9. Plot 2: right summary box
p2 <- ggplot(phylum_totals, aes(y = Phylum)) +
  geom_text(aes(x = 1, label = right_box_text), size = 3.8, fontface = "bold", color = "#2b2b2b") +
  labs(title = "Integrons\n(% of total)", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.background = element_rect(color = "#333333", fill = NA, linewidth = 1.0),
        plot.margin = margin(5, 10, 5, 10))

# 10. Combine + save
p3 <- ggplot(phylum_totals, aes(y = Phylum)) +
  geom_text(aes(x = 1, label = samp_text), size = 3.8, fontface = "bold", color = "#2b2b2b") +
  labs(title = "Samples with\nintegrons", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.background = element_rect(color = "#333333", fill = NA, linewidth = 1.0),
        plot.margin = margin(5, 10, 5, 10))

final_figure <- p1 + p2 + p3 + plot_layout(widths = c(4.0, 1.25, 1.25))
fig_h <- max(4.5, 1.5 + 0.75 * length(lv))
if (interactive()) print(final_figure)

ggsave(paste0(output_dir, "Fig2D_Integron_Phylum_Distribution.pdf"), final_figure,
       width = 14, height = fig_h, device = "pdf")
ggsave(paste0(output_dir, "Fig2D_Integron_Phylum_Distribution_600DPI.png"), final_figure,
       width = 14, height = fig_h, dpi = 600, bg = "white")
cat("\nFigure saved in", output_dir, "\n")
