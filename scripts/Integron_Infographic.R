# ==============================================================================
# INTEGRON HOST INFOGRAPHIC (R / ggplot2 + ggforce)
# KPI cards + host-phylum donut + top host genera (lollipop)
# Input : Integron_Phylum_Association.xlsx  (from Integron_Phylum_Pipeline.sh)
# Output: Fig_Integron_Hosts_Infographic.png (300 DPI) + .pdf
# The whole figure is one 14 x 8 inch canvas; positions below are in inches.
# ==============================================================================

# 1. Packages ------------------------------------------------------------------
pkgs <- c("readxl", "ggplot2", "dplyr", "ggforce")
new  <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) install.packages(new)
suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(dplyr); library(ggforce)
})

# 2. Paths ---------------------------------------------------------------------
in_dir  <- "Integron_Phylum_results/"
out_dir <- in_dir
EXC     <- c("Chordata", "Uroviricota")   # non-prokaryotic Kraken2 hits (host DNA, phage)

# 3. Data ----------------------------------------------------------------------
xlsx <- paste0(in_dir, "Integron_Phylum_Association.xlsx")
df <- read_excel(xlsx, sheet = "Integron_Phylum_Association")
ps <- read_excel(xlsx, sheet = "Per_sample")

n_all   <- nrow(df)
n_s_pos <- sum(ps$total > 0); n_s <- nrow(ps)
cls     <- df %>% filter(Phylum != "Unclassified", !Phylum %in% EXC)
known   <- !grepl("^Unclassified", cls$Organism_Name)
n_gen   <- n_distinct(cls$Organism_Name[known])
n_ph    <- n_distinct(cls$Phylum)

# 4. Style ---------------------------------------------------------------------
BG <- "#FFFBF5"; INK <- "#2B2D42"; MUTED <- "#8D99AE"
C_MAIN <- "#F28C6B"; C_BAC <- "#7BC8A4"; C_OTH <- "#A79AD8"; C_UNC <- "#E3E1DC"; C_BLUE <- "#6CB4EE"
fs <- function(points) points / .pt                       # font size in pt -> ggplot size
rect <- function(x, y, w, h, id) data.frame(x = c(x, x + w, x + w, x), y = c(y, y, y + h, y + h), id = id)

# 5. Donut data ----------------------------------------------------------------
n_main <- sum(df$Phylum == "Pseudomonadota")
n_bac  <- sum(df$Phylum == "Bacillota")
n_unc  <- sum(df$Phylum == "Unclassified")
n_oth  <- n_all - n_main - n_bac - n_unc
n_oth_ph <- n_distinct(df$Phylum[!df$Phylum %in% c("Pseudomonadota", "Bacillota", "Unclassified", EXC)])

don <- data.frame(lab = c("Pseudomonadota", "Bacillota", sprintf("Other phyla (%d)", n_oth_ph), "Unclassified"),
                  n   = c(n_main, n_bac, n_oth, n_unc),
                  col = c(C_MAIN, C_BAC, C_OTH, C_UNC))
don$end   <- cumsum(don$n) / n_all * 2 * pi          # 0 = 12 o'clock, clockwise
don$start <- c(0, head(don$end, -1))
cx <- 2.38; cy <- 2.96; R <- 1.9; R0 <- R * 0.66

leg <- don %>% mutate(i = row_number() - 1, yy = 3.95 - i * 0.78,
                      sub = sprintf("%s integrons  \u00b7  %.1f%%", prettyNum(n, big.mark = ","), 100 * n / n_all))
leg_sq <- do.call(rbind, lapply(seq_len(nrow(leg)), function(i)
  cbind(rect(4.75, leg$yy[i] - 0.11, 0.24, 0.24, i), col = leg$col[i])))

# 6. Top 12 genera -------------------------------------------------------------
gen <- cls[known, ] %>%
  group_by(Organism_Name) %>%
  summarise(total = n(), complete = sum(integron_type == "complete"),
            phy = first(Phylum), .groups = "drop") %>%
  arrange(desc(total), Organism_Name) %>% slice_head(n = 12)
ng <- nrow(gen)
x_l <- 9.31; x_r <- 13.58; y_b <- 0.80; y_t <- 5.08          # lollipop panel (inches)
xmax <- max(gen$total) * 1.5
xx <- function(v) x_l + v / xmax * (x_r - x_l)
gen$y <- y_b + ((ng - 1 - (seq_len(ng) - 1)) + 0.7) / ng * (y_t - y_b)
ticks <- pretty(c(0, xmax)); ticks <- ticks[ticks <= xmax]
all_main <- all(gen$phy == "Pseudomonadota")

# 7. KPI cards -----------------------------------------------------------------
cards <- data.frame(v = c(format(n_all, big.mark = ","), paste0(n_s_pos, "/", n_s), n_ph, n_gen),
                    l = c("integrons detected", "samples carry integrons", "host phyla", "host genera"),
                    col = c(C_MAIN, C_BLUE, C_BAC, C_OTH))
cw <- 2.75; ch <- 0.85; gap <- 0.28; y0 <- 6.1
x0 <- (14 - (4 * cw + 3 * gap)) / 2
cards$x <- x0 + (seq_len(4) - 1) * (cw + gap)
card_poly <- do.call(rbind, lapply(1:4, function(i) rect(cards$x[i], y0, cw, ch, i)))
bar_poly  <- do.call(rbind, lapply(1:4, function(i)
  cbind(rect(cards$x[i] + 0.16, y0 + 0.18, 0.07, ch - 0.36, i), col = cards$col[i])))

# 8. Build the canvas ----------------------------------------------------------
p <- ggplot() +
  coord_fixed(xlim = c(0, 14), ylim = c(0, 8), expand = FALSE) +
  scale_fill_identity() +
  theme_void() +
  theme(plot.background = element_rect(fill = BG, colour = NA), plot.margin = margin(0, 0, 0, 0)) +

  # title
  annotate("text", x = 7, y = 7.62, label = "Who carries integrons in the metagenome?",
           size = fs(24), fontface = "bold", colour = INK) +
  annotate("text", x = 7, y = 7.28,
           label = "Each integron linked to the taxonomy of its own contig  \u00b7  IntegronFinder 2 + Kraken2",
           size = fs(12.5), colour = MUTED) +

  # KPI cards
  geom_shape(data = card_poly, aes(x, y, group = id), radius = unit(3.5, "mm"),
             fill = "white", colour = "#EDE4D8", linewidth = 0.5) +
  geom_shape(data = bar_poly, aes(x, y, group = id, fill = col), radius = unit(0.8, "mm")) +
  geom_text(data = cards, aes(x + 0.4, y0 + 0.55, label = v), hjust = 0, size = fs(23),
            fontface = "bold", colour = INK) +
  geom_text(data = cards, aes(x + 0.4, y0 + 0.22, label = l), hjust = 0, size = fs(10.5), colour = MUTED) +

  # section titles
  annotate("text", x = 0.6, y = 5.45, label = "Host phylum of all integrons", hjust = 0, vjust = 0,
           size = fs(14.5), fontface = "bold", colour = INK) +
  annotate("text", x = 7.6, y = 5.45, label = "Top 12 host genera", hjust = 0, vjust = 0,
           size = fs(14.5), fontface = "bold", colour = INK) +
  { if (all_main) annotate("text", x = 7.6, y = 5.17, label = "all belong to Pseudomonadota", hjust = 0,
                           vjust = 0, size = fs(10.5), fontface = "bold", colour = C_MAIN) } +

  # donut
  geom_arc_bar(data = don, aes(x0 = cx, y0 = cy, r0 = R0, r = R, start = start, end = end, fill = col),
               colour = BG, linewidth = 1.2) +
  annotate("text", x = cx, y = cy + 0.19, label = sprintf("%.1f%%", 100 * n_main / n_all),
           size = fs(28), fontface = "bold", colour = C_MAIN) +
  annotate("text", x = cx, y = cy - 0.30, label = "Pseudomonadota",
           size = fs(11.5), fontface = "bold", colour = INK) +

  # donut legend
  geom_shape(data = leg_sq, aes(x, y, group = id, fill = col), radius = unit(1.5, "mm")) +
  geom_text(data = leg, aes(5.12, yy + 0.03, label = lab), hjust = 0, size = fs(11.5),
            fontface = "bold", colour = INK) +
  geom_text(data = leg, aes(5.12, yy - 0.25, label = sub), hjust = 0, size = fs(10), colour = MUTED) +

  # lollipop
  geom_segment(data = gen, aes(x = xx(0), xend = xx(total), y = y, yend = y),
               colour = "#F6C9B8", linewidth = 1.5, lineend = "round") +
  geom_point(data = gen, aes(xx(total), y), shape = 21, size = 5.6, fill = C_MAIN,
             colour = "white", stroke = 1.1) +
  geom_text(data = gen, aes(xx(total) + 0.26, y, label = total), hjust = 0, size = fs(11),
            fontface = "bold", colour = INK) +
  geom_text(data = filter(gen, complete > 0), aes(xx(total) + 0.68, y, label = paste(complete, "complete")),
            hjust = 0, size = fs(9), colour = "#3182bd") +
  geom_text(data = gen, aes(x_l - 0.12, y, label = Organism_Name), hjust = 1, size = fs(11.5),
            fontface = "italic", colour = INK) +
  annotate("segment", x = x_l, xend = x_r, y = y_b, yend = y_b, colour = "#D8D2C8", linewidth = 0.4) +
  annotate("segment", x = xx(ticks), xend = xx(ticks), y = y_b, yend = y_b - 0.05, colour = "#D8D2C8") +
  annotate("text", x = xx(ticks), y = y_b - 0.17, label = ticks, size = fs(10), colour = MUTED) +
  annotate("text", x = (x_l + x_r) / 2, y = y_b - 0.47, label = "Number of integrons",
           size = fs(11), colour = MUTED) +

  # footnote
  annotate("text", x = 0.6, y = 0.22, hjust = 0, size = fs(8.5), colour = MUTED,
           label = sprintf(paste0("Unclassified = contigs too short for k-mer classification (median ~0.8 kb).  ",
                                  "'Other phyla' includes %d non-prokaryotic hits.  Data unpublished; samples anonymised."),
                           sum(df$Phylum %in% EXC)))

print(p)
ggsave(paste0(out_dir, "Fig_Integron_Hosts_Infographic.png"), p, width = 14, height = 8, dpi = 300, bg = BG)
ggsave(paste0(out_dir, "Fig_Integron_Hosts_Infographic.pdf"), p, width = 14, height = 8, bg = BG)
cat("\nSaved:", paste0(out_dir, "Fig_Integron_Hosts_Infographic.png"), "\n")
