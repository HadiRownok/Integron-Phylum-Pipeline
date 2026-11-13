#!/usr/bin/env bash
# ==============================================================================
# Integron_Phylum_Pipeline.sh  (v3)
# Detect integrons in assembled metagenomic contigs and link every integron to
# the taxonomy of the contig that carries it (phylum / genus / species).
#
#   contigs --IntegronFinder 2--> integrons --Kraken2 on integron contigs-->
#   phylum + genus per integron --> tables (TSV + Excel) --> figure (R)
#   optional (-a): Kraken2 on ALL contigs --> DNA per phylum --> integron
#   enrichment per phylum (observed vs expected, binomial test) --> figure (R)
#
# Idea adapted from: Néron B. et al. (2022) IntegronFinder 2.0. Microorganisms
# 10:700, Fig. 2D (integron distribution by NCBI phylum). Pipeline code is ours.
#
# Usage:
#   conda activate <env with integron_finder + kraken2>
#   Integron_Phylum_Pipeline.sh -k /path/to/kraken2_db [options] *_contigs.fasta
#
# Options:
#   -k DIR    Kraken2 database folder (or set env KRAKEN2_DB)        [required*]
#   -e ENV    conda env that has kraken2, if not the active one
#   -t FILE   use a manual sample taxonomy table instead of Kraken2
#             (sample<TAB>Phylum<TAB>Organism_Name)                   [*alternative]
#   -o DIR    output folder                    (default: ./Integron_Phylum_results)
#   -r DIR    existing IntegronFinder results  (default: OUT/IntegronFinder)
#   -c N      CPUs per IntegronFinder run / Kraken2 threads   (default: 8)
#   -j N      samples run in parallel                         (default: 4)
#   -a        also run the integron-phylum enrichment analysis (Kraken2 on all contigs)
#   -h        help
#
# Sample name = file name without extension and without "_contigs"
#   sample1_contigs.fasta -> sample1
# Re-running is safe: finished IntegronFinder / Kraken2 steps are skipped.
# ==============================================================================
set -uo pipefail

KDB="${KRAKEN2_DB:-}"; KENV=""; TAXA=""; OUT="Integron_Phylum_results"
IF_DIR=""; CPU=8; JOBS=4; ENRICH=0

usage() { sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts ":k:e:t:o:r:c:j:ah" opt; do
    case $opt in
        k) KDB=$OPTARG ;;   e) KENV=$OPTARG ;;  t) TAXA=$OPTARG ;;
        o) OUT=$OPTARG ;;   r) IF_DIR=$OPTARG ;;
        c) CPU=$OPTARG ;;   j) JOBS=$OPTARG ;;  a) ENRICH=1 ;;
        h) usage 0 ;;
        :) echo "ERROR: -$OPTARG needs a value"; exit 1 ;;
        *) echo "ERROR: unknown option -$OPTARG"; usage 1 ;;
    esac
done
shift $((OPTIND - 1))

# ---------------------------------------------------------------- checks ----
die() { echo "ERROR: $*" >&2; exit 1; }
[[ $# -eq 0 ]] && die "give contig FASTA files, e.g.  $(basename "$0") -k DB *_contigs.fasta"
command -v integron_finder >/dev/null || die "integron_finder not found -> conda activate <env>"

if [[ -n $KENV ]]; then KRAKEN=(conda run --no-capture-output -n "$KENV" kraken2)
else KRAKEN=(kraken2); fi

if [[ -n $TAXA ]]; then
    [[ -s $TAXA ]] || die "taxonomy file '$TAXA' not found"
    MODE=manual
else
    [[ -n $KDB ]] || die "give a Kraken2 database with -k (or a manual table with -t)"
    [[ -f $KDB/hash.k2d ]] || die "'$KDB' is not a Kraken2 database (hash.k2d missing)"
    "${KRAKEN[@]}" --version >/dev/null 2>&1 || die "kraken2 not found (activate its env or use -e ENV)"
    MODE=kraken2
fi

IF_DIR="${IF_DIR:-$OUT/IntegronFinder}"
K_DIR="$OUT/Kraken2"
LOG="$OUT/pipeline.log"
mkdir -p "$OUT" "$IF_DIR"
echo "=== $(date) | $# files | taxonomy: $MODE | CPU=$CPU | parallel=$JOBS ===" | tee -a "$LOG"

sample_of() { local b; b=$(basename "$1"); b=${b%.*}; echo "${b%_contigs}"; }
base_of()   { local b; b=$(basename "$1"); echo "${b%.*}"; }
is_done()   { ls "$IF_DIR/$1/Results_Integron_Finder_$1/"*.summary >/dev/null 2>&1; }

# ==============================================================================
# STEP 1  IntegronFinder 2 on every contig file
# ==============================================================================
run_if() {
    local fa=$1 base sdir
    base=$(base_of "$fa"); sdir="$IF_DIR/$base"
    if is_done "$base"; then echo "[skip] IntegronFinder $base" | tee -a "$LOG"; return; fi
    mkdir -p "$sdir"
    echo "[run ] IntegronFinder $base ($(date +%H:%M))" | tee -a "$LOG"
    if integron_finder "$fa" --cpu "$CPU" --func-annot --outdir "$sdir" \
                       --pdf --gbk --mute > "$sdir/run.log" 2>&1; then
        echo "[ ok ] $base" | tee -a "$LOG"
    else
        echo "[FAIL] $base (see $sdir/run.log)" | tee -a "$LOG"
    fi
}

for fa in "$@"; do
    [[ -s $fa ]] || { echo "[warn] missing/empty: $fa" | tee -a "$LOG"; continue; }
    run_if "$fa" &
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done
done
wait

# ==============================================================================
# STEP 2  One row per integron from *.integrons
#   cols: 1 ID_integron 2 ID_replicon 8 type_elt 9 annotation 11 type
#   #attC = annotation "attC" ; #protein = type_elt "protein" (intI included)
# ==============================================================================
RAW="$OUT/Integron_per_contig_raw.tsv"
printf "sample\treplicon_name\tintegron_in_replicon_id\tintegron_type\t#attC\t#protein\n" > "$RAW"
for fa in "$@"; do
    base=$(base_of "$fa"); s=$(sample_of "$fa")
    for f in "$IF_DIR/$base/Results_Integron_Finder_$base/"*.integrons; do
        [[ -e $f ]] || continue
        awk -F'\t' -v s="$s" '
            /^#/ || /^ID_integron/ || NF < 11 { next }
            { k = $2 SUBSEP $1
              if (!(k in t)) { o[++n] = k; t[k] = $11 }
              if ($9 == "attC")    a[k]++
              if ($8 == "protein") p[k]++ }
            END { for (i = 1; i <= n; i++) { split(o[i], x, SUBSEP)
                  printf "%s\t%s\t%s\t%s\t%d\t%d\n", s, x[1], x[2], t[o[i]], a[o[i]], p[o[i]] } }' "$f"
    done
done | sort -t$'\t' -k1,1V -k2,2 -k3,3 >> "$RAW"
N_INT=$(($(wc -l < "$RAW") - 1))
echo "[done] $N_INT integrons found" | tee -a "$LOG"
(( N_INT > 0 )) || { echo "No integrons in any sample. Stopping." | tee -a "$LOG"; exit 0; }

# ==============================================================================
# STEP 3  Taxonomy for every integron-carrying contig
# ==============================================================================
CTAX="$OUT/contig_taxonomy.tsv"     # sample  contig  Phylum  Genus  Species  taxid  kraken_name

if [[ $MODE == kraken2 ]]; then
    mkdir -p "$K_DIR"
    KFA="$K_DIR/integron_contigs.fasta"
    KOUT="$K_DIR/integron_contigs.kraken"
    KREP="$K_DIR/integron_contigs.kreport"

    # 3a. extract only the contigs that carry an integron; header = sample|contig
    : > "$KFA.tmp"
    for fa in "$@"; do
        s=$(sample_of "$fa")
        awk -F'\t' -v s="$s" 'NR > 1 && $1 == s { print $2 }' "$RAW" | sort -u > "$K_DIR/.ids"
        [[ -s "$K_DIR/.ids" ]] || continue
        awk -v s="$s" 'NR == FNR { want[$1]; next }
            /^>/ { id = substr($1, 2); keep = (id in want); if (keep) print ">" s "|" id; next }
            keep' "$K_DIR/.ids" "$fa" >> "$KFA.tmp"
    done
    rm -f "$K_DIR/.ids"
    # keep the old file (and its timestamp) if nothing changed -> Kraken2 is skipped
    if cmp -s "$KFA.tmp" "$KFA"; then rm "$KFA.tmp"; else mv "$KFA.tmp" "$KFA"; fi
    echo "[done] $(grep -c '^>' "$KFA") integron contigs extracted -> $KFA" | tee -a "$LOG"

    # 3b. Kraken2 (one run for all samples, so the database loads only once)
    if [[ -s $KREP && $KREP -nt $KFA ]]; then
        echo "[skip] Kraken2 (report is newer than contigs)" | tee -a "$LOG"
    else
        echo "[run ] Kraken2 on integron contigs ($(date +%H:%M))" | tee -a "$LOG"
        "${KRAKEN[@]}" --db "$KDB" --threads "$CPU" --use-names \
            --report "$KREP" --output "$KOUT" "$KFA" 2>> "$LOG" \
            || die "Kraken2 failed, see $LOG"
    fi

    # 3c. lineage of every taxid from the report tree (indentation = depth),
    #     then phylum / genus / species for every contig
    awk -F'\t' -v OFS='\t' '
        NR == FNR {                                  # kraken2 report
            rank = $(NF-2); tid = $(NF-1); nm = $NF
            raw = nm; sub(/^ +/, "", nm); d = (length(raw) - length(nm)) / 2 + 1
            P[d] = (rank == "P") ? nm : P[d-1]
            G[d] = (rank == "G") ? nm : (rank ~ /^[DKPCOF]/ ? "" : G[d-1])
            S[d] = (rank == "S") ? nm : (rank ~ /^S[0-9]/ ? S[d-1] : "")
            ph[tid] = P[d]; ge[tid] = G[d]; sp[tid] = S[d]
            next }
        {                                            # kraken2 per-contig output
            split($2, id, "|")
            tid = $3; sub(/.*\(taxid /, "", tid); sub(/\).*/, "", tid)
            name = $3; sub(/ \(taxid [0-9]+\)$/, "", name)
            if ($1 == "U") { p = "Unclassified"; g = "Unclassified"; s = "" }
            else { p = (ph[tid] != "" ? ph[tid] : "Unclassified")
                   g = (ge[tid] != "" ? ge[tid] : "Unclassified_" p); s = sp[tid] }
            print id[1], id[2], p, g, (s == "" ? "NA" : s), tid, name }' \
        "$KREP" "$KOUT" > "$CTAX"
else
    # manual: same genus/phylum for every contig of a sample
    awk -F'\t' -v OFS='\t' 'NR == FNR { if (FNR > 1) t[$1] = $2 OFS $3; next }
        FNR > 1 { split(($1 in t ? t[$1] : "NA" OFS "NA"), x, OFS)
                  print $1, $2, x[1], x[2], "NA", "NA", "manual" }' "$TAXA" "$RAW" > "$CTAX"
fi

# ==============================================================================
# STEP 4  Main table = integron rows + contig taxonomy
# ==============================================================================
MAIN="$OUT/Integron_Phylum_Association.tsv"
awk -F'\t' -v OFS='\t' '
    NR == FNR { tx[$1 SUBSEP $2] = $3 OFS $4 OFS $5 OFS $6 OFS $7; next }
    FNR == 1  { print $0, "Phylum", "Organism_Name", "Species", "taxid", "kraken_name"; next }
    { k = $1 SUBSEP $2; print $0, (k in tx ? tx[k] : "Unclassified" OFS "Unclassified" OFS "NA" OFS "0" OFS "NA") }' \
    "$CTAX" "$RAW" > "$MAIN"

# ==============================================================================
# STEP 5  Summaries
# ==============================================================================
PER_SAMPLE="$OUT/Integron_counts_per_sample.tsv"
{
  printf "sample\tstatus\tcomplete\tIn0\tCALIN\ttotal\n"
  for fa in "$@"; do
      base=$(base_of "$fa"); s=$(sample_of "$fa")
      is_done "$base" && st=done || st=failed
      awk -F'\t' -v s="$s" -v st="$st" 'FNR > 1 && $1 == s { c[$4]++ }
          END { printf "%s\t%s\t%d\t%d\t%d\t%d\n", s, st, c["complete"], c["In0"], c["CALIN"], c["complete"]+c["In0"]+c["CALIN"] }' "$RAW"
  done | sort -V
} > "$PER_SAMPLE"

summ() {   # $1 = column of MAIN to group by, $2 = header name
    awk -F'\t' -v col="$1" -v h="$2" '
        FNR > 1 { g = $col; n[g]++; c[g, $4]++; N++ }
        END { printf "%s\tcomplete\tIn0\tCALIN\ttotal\tpercent\n", h
              for (g in n) printf "%s\t%d\t%d\t%d\t%d\t%.1f\n", g, c[g,"complete"], c[g,"In0"], c[g,"CALIN"], n[g], 100*n[g]/N }' "$MAIN" \
    | { read -r hd; echo "$hd"; sort -t$'\t' -k5,5nr; }
}
PER_PHYLUM="$OUT/Integron_counts_per_phylum.tsv"; summ 7 Phylum        > "$PER_PHYLUM"
PER_GENUS="$OUT/Integron_counts_per_genus.tsv";   summ 8 Organism_Name > "$PER_GENUS"

# ==============================================================================
# STEP 6  Excel workbook (needs python3 + pandas + openpyxl)
# ==============================================================================
XLSX="$OUT/Integron_Phylum_Association.xlsx"
python3 - "$XLSX" "$MAIN" "$PER_SAMPLE" "$PER_PHYLUM" "$PER_GENUS" <<'PY' 2>/dev/null \
    && echo "[done] Excel -> $XLSX" | tee -a "$LOG" \
    || echo "[warn] Excel not written (needs pandas + openpyxl); TSV files are complete" | tee -a "$LOG"
import sys, pandas as pd
x, *files = sys.argv[1:]
names = ["Integron_Phylum_Association", "Per_sample", "Per_phylum", "Per_genus"]
with pd.ExcelWriter(x) as w:
    for f, n in zip(files, names):
        pd.read_csv(f, sep="\t", keep_default_na=False).to_excel(w, sheet_name=n, index=False)
PY

# ==============================================================================
# STEP 7  Figure (R: ggplot2, dplyr, tidyr, patchwork). R scripts saved in OUT/R/
# ==============================================================================
mkdir -p "$OUT/R"
RFIG="$OUT/R/Integron_Phylum_Figure.R"
REN="$OUT/R/Integron_Phylum_Enrichment.R"
cat > "$RFIG" <<'RS'
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
  "Integron_Phylum_Association.tsv"
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
RS
cat > "$REN" <<'RS'
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
RS

if command -v Rscript >/dev/null && Rscript "$RFIG" "$MAIN" "$OUT" "$#" >> "$LOG" 2>&1; then
    echo "[done] Figure -> $OUT/Fig2D_Integron_Phylum_Distribution_600DPI.png" | tee -a "$LOG"
else
    echo "[warn] Figure not made here (needs Rscript + ggplot2 dplyr tidyr patchwork; see $LOG)." | tee -a "$LOG"
    echo "       In RStudio: source $RFIG after setting file_path to $MAIN" | tee -a "$LOG"
fi

# ==============================================================================
# STEP 8  (optional, -a) Integron-phylum enrichment: observed vs expected
#   Kraken2 on ALL contigs -> bp of DNA per phylum per sample
# ==============================================================================
if (( ENRICH )); then
  if [[ $MODE != kraken2 ]]; then
    echo "[warn] -a needs Kraken2 (-k); enrichment skipped" | tee -a "$LOG"
  else
    KA="$OUT/Kraken2_all_contigs"; mkdir -p "$KA"
    for fa in "$@"; do
        s=$(sample_of "$fa")
        if [[ -s $KA/$s.taxid_bp.tsv && -s $KA/$s.kreport ]]; then
            echo "[skip] Kraken2 all contigs $s" | tee -a "$LOG"; continue
        fi
        echo "[run ] Kraken2 all contigs $s ($(date +%H:%M))" | tee -a "$LOG"
        # write per-contig output to a file (stdout does not pass reliably through "conda run")
        "${KRAKEN[@]}" --db "$KDB" --threads "$CPU" --memory-mapping \
            --report "$KA/$s.kreport" --output "$KA/$s.kraken" "$fa" 2>> "$LOG" \
        && awk -F'\t' '{ bp[$3] += $4; n[$3]++ } END { for (t in bp) print t "\t" bp[t] "\t" n[t] }' \
            "$KA/$s.kraken" > "$KA/$s.taxid_bp.tmp" \
        && mv "$KA/$s.taxid_bp.tmp" "$KA/$s.taxid_bp.tsv" && rm -f "$KA/$s.kraken" \
        || { echo "[FAIL] Kraken2 all contigs $s" | tee -a "$LOG"; rm -f "$KA/$s.taxid_bp.tmp"; }
    done

    BP="$OUT/Phylum_bp_per_sample.tsv"
    printf "sample\tPhylum\tbp\tn_contigs\n" > "$BP"
    for fa in "$@"; do
        s=$(sample_of "$fa")
        [[ -s $KA/$s.taxid_bp.tsv ]] || continue
        awk -F'\t' -v OFS='\t' -v s="$s" '
            NR == FNR {
                rank = $(NF-2); tid = $(NF-1); nm = $NF
                raw = nm; sub(/^ +/, "", nm); d = (length(raw) - length(nm)) / 2 + 1
                P[d] = (rank == "P") ? nm : P[d-1]; ph[tid] = P[d]; next }
            { p = ($1 == 0 || ph[$1] == "") ? "Unclassified" : ph[$1]; bp[p] += $2; n[p] += $3 }
            END { for (p in bp) print s, p, bp[p], n[p] }' "$KA/$s.kreport" "$KA/$s.taxid_bp.tsv"
    done | sort -t$'\t' -k1,1V -k3,3nr >> "$BP"

    EIN="$OUT/Integron_Phylum_Enrichment_input.tsv"
    awk -F'\t' -v OFS='\t' '
        NR == FNR { if (FNR > 1) { k = $1 SUBSEP $7; ni[k]++; nc[k] += ($4 == "complete") }; next }
        FNR == 1  { print "sample", "Phylum", "bp", "n_contigs", "integrons", "complete"; next }
        { k = $1 SUBSEP $2; print $1, $2, $3, $4, (k in ni ? ni[k] : 0), (k in nc ? nc[k] : 0) }' \
        "$MAIN" "$BP" > "$EIN"
    echo "[done] enrichment input -> $EIN" | tee -a "$LOG"

    if command -v Rscript >/dev/null && Rscript "$REN" "$OUT" >> "$LOG" 2>&1; then
        echo "[done] Enrichment -> $OUT/Integron_Phylum_Enrichment_results.csv + Fig_Integron_Phylum_Enrichment_600DPI.png" | tee -a "$LOG"
    else
        echo "[warn] Enrichment statistics/figure not made here (needs R + ggplot2 dplyr patchwork ggrepel scales)." | tee -a "$LOG"
        echo "       Copy $EIN to a PC and run $REN in RStudio." | tee -a "$LOG"
    fi
  fi
fi

# ------------------------------------------------------------------ report --
echo "" | tee -a "$LOG"
echo "[summary] integrons: $N_INT | samples with integrons: $(awk -F'\t' 'NR>1 && $6>0' "$PER_SAMPLE" | wc -l)/$(($(wc -l < "$PER_SAMPLE") - 1))" | tee -a "$LOG"
awk -F'\t' 'NR>1 && $2=="failed" { print "[FAIL] " $1 ": no IntegronFinder result" }' "$PER_SAMPLE" | tee -a "$LOG"
{ command -v column >/dev/null && column -t -s$'\t' "$PER_PHYLUM" || cat "$PER_PHYLUM"; } | tee -a "$LOG"
echo "Outputs: $OUT/" | tee -a "$LOG"
