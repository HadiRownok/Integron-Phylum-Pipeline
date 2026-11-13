# ==============================================================================
# INTEGRON OVERVIEW FIGURES (multi-panel, publication / LinkedIn quality)
# Input : Integron_Phylum_Association.xlsx  (from Integron_Phylum_Pipeline.sh)
# Output: Fig_Integron_Overview.pdf / _600DPI.png  + every panel as its own PNG
#
#   A  Integrons per sample (complete / In0 / CALIN), all samples incl. zero
#   B  Sample x phylum heatmap (where the integrons come from)
#   C  Top 15 host genera (Kraken2, integron-carrying contigs)
#   D  Cassette array size: attC sites per integron (complete vs CALIN)
#   E  Contig length of classified vs unclassified integron contigs (Wilcoxon)
# ==============================================================================

# 1. Packages ------------------------------------------------------------------
pkgs <- c("readxl", "ggplot2", "dplyr", "tidyr", "patchwork", "scales")
new  <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) install.packages(new)
suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(dplyr)
  library(tidyr);  library(patchwork); library(scales)
})

# 2. Paths ---------------------------------------------------------------------
# Folder that contains Integron_Phylum_Association.xlsx (pipeline output)
in_dir  <- "Integron_Phylum_results/"
xlsx    <- paste0(in_dir, "Integron_Phylum_Association.xlsx")
out_dir <- in_dir
fig_title <- "Integrons in metagenomes: abundance, taxonomic origin and cassette structure"

# 3. Data ----------------------------------------------------------------------
df <- read_excel(xlsx, sheet = "Integron_Phylum_Association")
ps <- read_excel(xlsx, sheet = "Per_sample")

type_levels <- c("complete", "In0", "CALIN")                  # legend order
type_colors <- c(complete = "#3182bd", In0 = "#de2d26", CALIN = "#fec44f")
df$integron_type <- factor(df$integron_type, levels = rev(type_levels))  # CALIN stacks first

sample_order <- ps$sample[order(as.numeric(gsub("\\D", "", ps$sample)))]

# Anonymise sample IDs for figures (unpublished / confidential data)
anonymize <- TRUE
if (anonymize) {
  new_id <- setNames(paste0("S", seq_along(sample_order)), sample_order)
  df$sample <- unname(new_id[df$sample]); ps$sample <- unname(new_id[ps$sample])
  sample_order <- unname(new_id[sample_order])
}

# phyla < 1% of integrons -> "Other"
ph_n  <- df %>% count(Phylum) %>% mutate(pct = 100 * n / sum(n))
minor <- ph_n$Phylum[ph_n$pct < 1 & ph_n$Phylum != "Unclassified"]
df$Phylum_g <- ifelse(df$Phylum %in% minor, "Other", df$Phylum)
main_ph  <- ph_n %>% filter(!Phylum %in% c(minor, "Unclassified")) %>% arrange(n) %>% pull(Phylum)
ph_levels <- c("Unclassified", if (length(minor)) "Other", main_ph)  # bottom -> top

# contig length from SPAdes name: NODE_12_length_3456_cov_7.8
df$contig_len <- as.numeric(sub(".*_length_([0-9]+)_.*", "\\1", df$replicon_name))
df$tax_status <- ifelse(df$Phylum == "Unclassified", "Unclassified", "Classified")

# 4. Common theme --------------------------------------------------------------
theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, colour = "grey25"),
        axis.text     = element_text(colour = "black"),
        axis.title    = element_text(face = "bold"),
        legend.position = "top",
        legend.title  = element_text(face = "bold", size = 9.5),
        legend.text   = element_text(size = 9),
        plot.tag      = element_text(face = "bold", size = 15))

fill_type <- scale_fill_manual(values = type_colors, breaks = type_levels, name = "Integron type")

# 5. A: integrons per sample ---------------------------------------------------
ps_long <- ps %>%
  select(sample, complete, In0, CALIN) %>%
  pivot_longer(-sample, names_to = "integron_type", values_to = "n") %>%
  mutate(sample = factor(sample, levels = sample_order),
         integron_type = factor(integron_type, levels = rev(type_levels)))
zero <- ps %>% filter(total == 0) %>% mutate(sample = factor(sample, levels = sample_order))
n_pos <- sum(ps$total > 0)

pA <- ggplot(ps_long, aes(sample, n, fill = integron_type)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.72, colour = "black", linewidth = 0.2) +
  geom_text(data = zero, aes(sample, 0, label = "0"), inherit.aes = FALSE, vjust = -0.5, size = 3) +
  fill_type +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(title = "Integron load per metagenome",
       subtitle = sprintf("%d integrons; %d of %d samples carry at least one", nrow(df), n_pos, nrow(ps)),
       x = NULL, y = "Number of integrons") +
  theme_pub + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5))

# 6. B: sample x phylum heatmap ------------------------------------------------
hm <- df %>%
  count(sample, Phylum_g, name = "n") %>%
  complete(sample = sample_order, Phylum_g = ph_levels, fill = list(n = 0)) %>%
  mutate(sample = factor(sample, levels = sample_order),
         Phylum_g = factor(Phylum_g, levels = ph_levels))

pB <- ggplot(hm, aes(sample, Phylum_g, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(n > 0, n, ""), colour = n > 40), size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", trans = "log1p",
                      breaks = c(0, 1, 10, 100), name = "Integrons") +
  labs(title = "Taxonomic origin of integrons per sample",
       subtitle = "Kraken2 classification of integron-carrying contigs", x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5),
        axis.text.y = element_text(face = "bold"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.width = unit(1.2, "cm"))

# 7. C: top 15 host genera -----------------------------------------------------
gen <- df %>% filter(!grepl("^Unclassified", Organism_Name)) %>%
  count(Organism_Name, integron_type, name = "n")
top <- gen %>% group_by(Organism_Name) %>% summarise(t = sum(n), .groups = "drop") %>%
  slice_max(t, n = 15, with_ties = FALSE) %>% arrange(t)
gen <- gen %>% filter(Organism_Name %in% top$Organism_Name) %>%
  mutate(Organism_Name = factor(Organism_Name, levels = top$Organism_Name))
n_gen_total <- n_distinct(df$Organism_Name[!grepl("^Unclassified", df$Organism_Name)])

pC <- ggplot(gen, aes(n, Organism_Name, fill = integron_type)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7, colour = "black", linewidth = 0.2) +
  fill_type +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Top 15 host genera",
       subtitle = sprintf("of %d genera assigned at genus level", n_gen_total),
       x = "Number of integrons", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(face = "italic"))

# 8. D: cassette array size (attC per integron) -------------------------------
#    In0 has no attC by definition; CALIN requires >= 2 attC by definition,
#    so this is shown as a distribution, without a statistical test.
attc <- df %>% filter(integron_type %in% c("complete", "CALIN")) %>%
  count(integron_type, attC = `#attC`, name = "n") %>%
  mutate(integron_type = factor(integron_type, levels = c("complete", "CALIN")))
max_attc <- max(attc$attC)

pD <- ggplot(attc, aes(factor(attC, levels = 1:max_attc), n, fill = integron_type)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.75, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = type_colors[c("complete", "CALIN")], name = "Integron type") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 10, 50, 100, 500),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(drop = FALSE) +
  labs(title = "Cassette array size",
       subtitle = "attC sites per element (In0 carry none)",
       x = "Number of attC sites", y = "Integrons (log scale)") +
  theme_pub

# 9. E: contig length vs taxonomic assignment ---------------------------------
len <- df %>% distinct(sample, replicon_name, .keep_all = TRUE)
wt  <- suppressWarnings(wilcox.test(contig_len ~ tax_status, data = len))
med <- len %>% group_by(tax_status) %>% summarise(m = median(contig_len), k = n(), .groups = "drop")
p_txt <- if (wt$p.value < 1e-4) "p < 0.0001" else sprintf("p = %.2g", wt$p.value)

pE <- ggplot(len, aes(tax_status, contig_len, fill = tax_status)) +
  geom_violin(colour = "black", linewidth = 0.3, alpha = 0.85, scale = "width") +
  geom_boxplot(width = 0.14, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  geom_text(data = med, aes(tax_status, max(len$contig_len) * 1.6,
                            label = paste0("n = ", k, "\nmedian ", comma(m), " bp")),
            inherit.aes = FALSE, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = c(Classified = "#74a9cf", Unclassified = "#bdbdbd"), guide = "none") +
  scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "Why are some integrons unclassified?",
       subtitle = paste("Contig length, Wilcoxon rank-sum", p_txt),
       x = NULL, y = "Contig length (bp, log scale)") +
  theme_pub + theme(axis.text.x = element_text(face = "bold"))

# 10. Combine + save -----------------------------------------------------------
fig <- (pA | pB) / (pC | pD | pE) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    tag_levels = "A",
    title = fig_title,
    caption = "Integrons detected with IntegronFinder 2; host taxonomy of integron-carrying contigs by Kraken2 (Standard-16 DB).",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.caption = element_text(size = 9, colour = "grey30")))

ggsave(paste0(out_dir, "Fig_Integron_Overview.pdf"), fig, width = 18, height = 11, device = "pdf")
ggsave(paste0(out_dir, "Fig_Integron_Overview_600DPI.png"), fig, width = 18, height = 11, dpi = 600, bg = "white")

# each panel separately (LinkedIn carousel / slides)
panels <- list(A_per_sample = pA, B_heatmap = pB, C_top_genera = pC,
               D_attC_arrays = pD, E_contig_length = pE)
sizes  <- list(A_per_sample = c(9, 5), B_heatmap = c(9, 4.5), C_top_genera = c(6.5, 5.5),
               D_attC_arrays = c(6.5, 5), E_contig_length = c(5.5, 5.5))
for (nm in names(panels)) {
  ggsave(paste0(out_dir, "Panel_", nm, ".png"), panels[[nm]],
         width = sizes[[nm]][1], height = sizes[[nm]][2], dpi = 600, bg = "white")
}

print(fig)
cat("\nSaved in", out_dir, "\n")
cat(sprintf("Contig length: classified median %s bp vs unclassified %s bp, %s\n",
            comma(med$m[med$tax_status == "Classified"]),
            comma(med$m[med$tax_status == "Unclassified"]), p_txt))
