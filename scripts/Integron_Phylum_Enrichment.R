# ==============================================================================
# INTEGRON–PHYLUM ENRICHMENT: observed vs expected integrons per phylum
# Input : Integron_Phylum_Enrichment_input.tsv  (from Integron_Phylum_Enrichment.sh)
# Output: Integron_Phylum_Enrichment_results.csv
#         Fig_Integron_Phylum_Enrichment.pdf / _600DPI.png
#
# Expected integrons of a phylum = all classified integrons x that phylum's share
# of classified metagenomic DNA (bp). Tested with an exact binomial test,
# Benjamini–Hochberg correction across phyla.
# ==============================================================================

# Run in RStudio (edit in_dir below) or:  Rscript Integron_Phylum_Enrichment.R <folder>
args <- commandArgs(trailingOnly = TRUE)
pkgs <- c("ggplot2", "dplyr", "patchwork", "ggrepel", "scales")
new  <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) {
  if (!interactive()) stop("Missing R packages: ", paste(new, collapse = ", "))
  install.packages(new, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(ggrepel); library(scales)
})

# ---- settings -----------------------------------------------------------------
in_dir   <- if (length(args) >= 1) args[1] else "Integron_Phylum_results/"
in_dir   <- paste0(sub("/$", "", in_dir), "/")
out_dir  <- in_dir
EXC      <- c("Chordata", "Uroviricota")  # non-prokaryotic hits
min_share <- 0.005                         # show phyla with >= 0.5% of DNA ...
min_int   <- 2                             # ... or >= 2 integrons

# ---- data ---------------------------------------------------------------------
d <- read.delim(paste0(in_dir, "Integron_Phylum_Enrichment_input.tsv"), stringsAsFactors = FALSE) %>%
  filter(Phylum != "Unclassified", !Phylum %in% EXC)

pool <- d %>% group_by(Phylum) %>%
  summarise(bp = sum(bp), integrons = sum(integrons), complete = sum(complete),
            samples_present = sum(bp > 0), .groups = "drop")
N <- sum(pool$integrons); B <- sum(pool$bp)

pool <- pool %>% mutate(
  dna_share      = bp / B,
  integron_share = integrons / N,
  expected       = N * dna_share,
  per_Mb         = integrons / (bp / 1e6),
  enrichment     = integron_share / dna_share)

bt <- lapply(seq_len(nrow(pool)), function(i) binom.test(pool$integrons[i], N, pool$dna_share[i]))
pool$p_value  <- sapply(bt, `[[`, "p.value")
pool$ci_low   <- sapply(bt, function(x) x$conf.int[1]) / pool$dna_share
pool$ci_high  <- sapply(bt, function(x) x$conf.int[2]) / pool$dna_share
pool$q_value  <- p.adjust(pool$p_value, method = "BH")
pool$status   <- case_when(pool$q_value < 0.05 & pool$enrichment > 1 ~ "Enriched",
                           pool$q_value < 0.05 & pool$enrichment < 1 ~ "Depleted",
                           TRUE ~ "Not significant")
pool$stars    <- case_when(pool$q_value < 0.001 ~ "***", pool$q_value < 0.01 ~ "**",
                           pool$q_value < 0.05 ~ "*", TRUE ~ "")

write.csv(pool %>% arrange(desc(enrichment)), paste0(out_dir, "Integron_Phylum_Enrichment_results.csv"), row.names = FALSE)

# ---- plotting data ------------------------------------------------------------
show <- pool %>% filter(dna_share >= min_share | integrons >= min_int)
lim  <- 6                                                  # log2 axis limit
show <- show %>% mutate(
  l2   = pmin(pmax(log2((integrons + 0.5) / (expected + 0.5)), -lim), lim),   # 0.5 pseudocount
  lo   = pmax(log2(pmax(ci_low, 2^-lim)), -lim),
  hi   = pmin(log2(pmax(ci_high, 2^-lim)), lim),
  lab  = sprintf("%d obs / %.0f exp %s", integrons, expected, stars),
  Phylum = factor(Phylum, levels = Phylum[order(l2)]))
cols <- c(Enriched = "#E4572E", Depleted = "#2E86AB", `Not significant` = "#B8B8B8")

theme_pub <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        axis.text = element_text(colour = "black"), axis.title = element_text(face = "bold"),
        legend.position = "top", legend.title = element_blank(),
        plot.tag = element_text(face = "bold", size = 15))

# A. observed / expected (forest plot)
pA <- ggplot(show, aes(l2, Phylum, colour = status)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.6) +
  geom_point(aes(size = dna_share * 100)) +
  geom_text(aes(x = lim + 0.3, label = lab), hjust = 0, size = 3.2, colour = "grey15") +
  scale_colour_manual(values = cols) +
  scale_size_area(max_size = 8, name = "% of DNA", guide = guide_legend(order = 2)) +
  scale_x_continuous(limits = c(-lim, lim + 4.2), breaks = seq(-lim, lim, 2),
                     labels = function(x) ifelse(x > lim, "", x)) +
  labs(title = "Which phyla carry more integrons than expected?",
       subtitle = "Observed vs expected from each phylum's share of metagenomic DNA\nexact binomial test, BH-adjusted (* q<0.05, ** q<0.01, *** q<0.001)",
       x = expression(bold(log[2]~"(observed / expected integrons)")), y = NULL) +
  theme_pub + theme(axis.text.y = element_text(face = "bold"))

# B. DNA share vs integron share
floor_v <- 0.5 / N                                        # where zero-integron phyla are drawn
show <- show %>% mutate(int_plot = ifelse(integrons == 0, floor_v, integron_share))
rng  <- range(c(show$dna_share, show$int_plot))
pB <- ggplot(show, aes(dna_share, int_plot, colour = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = rng[1] * 1.5, y = rng[2] * 0.8, label = "more integrons\nthan expected",
           hjust = 0, vjust = 1, size = 3.3, colour = "#E4572E", fontface = "italic") +
  annotate("text", x = rng[2] * 0.8, y = rng[1] * 1.5, label = "fewer integrons\nthan expected",
           hjust = 1, vjust = 0, size = 3.3, colour = "#2E86AB", fontface = "italic") +
  geom_point(size = 3.4) +
  geom_text_repel(aes(label = Phylum), size = 3.2, colour = "grey15", max.overlaps = 30,
                  min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_log10(labels = label_percent(accuracy = 0.1)) +
  scale_y_log10(labels = function(x) ifelse(x <= floor_v * 1.01, "0", label_percent(accuracy = 0.1)(x))) +
  coord_equal() +
  labs(title = "Share of DNA vs share of integrons",
       subtitle = "Dashed line: integrons proportional to abundance",
       x = "Share of classified metagenomic DNA", y = "Share of classified integrons") +
  theme_pub

fig <- pA + pB + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "A")
h <- max(6, 2 + 0.38 * nrow(show))
ggsave(paste0(out_dir, "Fig_Integron_Phylum_Enrichment.pdf"), fig, width = 16, height = h)
ggsave(paste0(out_dir, "Fig_Integron_Phylum_Enrichment_600DPI.png"), fig, width = 16, height = h, dpi = 600, bg = "white")
if (interactive()) print(fig)
print(pool %>% arrange(q_value) %>% select(Phylum, integrons, expected, enrichment, q_value, status), n = 30)
