# ============================================================
# Volcano plot: EPI300_pCC1_ATF1 (noAra vs iAAl)
# Styled to match Met_EPPCC1_EPFOS_iAA_Volcano.png
# ============================================================

library(tidyverse)   # dplyr, ggplot2, readr
library(ggrepel)     # non-overlapping labels with leader lines

# ---- 1. Load data -------------------------------------------------
# Expects a CSV with columns: Name, EPI300_pCC1_ATF1_noAra, EPI300_pCC1_ATF1_iAAl
df <- read_csv("EPI300_pCC1_ATF1_noAraiAAlOD.csv", show_col_types = FALSE)

# ---- 2. Fold change -------------------------------------------------
# Many Names are detected in only one condition (blank cell in the CSV).
# Dropping them (or letting NA propagate through sum()/log2()) would remove
# points from the plot, so blanks are treated as 0 counts, and a small
# pseudocount is added ONLY inside the fold-change ratio so log2(0) and
# division by zero never occur. Every row in the file stays on the plot.
pseudo <- 0.5

df <- df %>%
  mutate(
    EPI300_pCC1_ATF1_noAra = replace_na(EPI300_pCC1_ATF1_noAra, 0),
    EPI300_pCC1_ATF1_iAAl  = replace_na(EPI300_pCC1_ATF1_iAAl,   0)
  )

total_noAra <- sum(df$EPI300_pCC1_ATF1_noAra)
total_iAAl   <- sum(df$EPI300_pCC1_ATF1_iAAl)

df <- df %>%
  mutate(
    prop_noAra = EPI300_pCC1_ATF1_noAra / total_noAra,
    prop_iAAl   = EPI300_pCC1_ATF1_iAAl  / total_iAAl,
    log2FC     = log2(((EPI300_pCC1_ATF1_iAAl  + pseudo) / total_iAAl) /
                      ((EPI300_pCC1_ATF1_noAra + pseudo) / total_noAra))
  )

# ---- 3. Significance ------------------------------------------------
# NOTE: this dataset has ONE measurement per Name per condition (no
# biological replicates), so a textbook p-value isn't really available.
# As a stand-in, this uses a two-proportion z-test of each Name's
# share of the total signal in Ara vs noAra. Treat these p-values as a
# rough, exploratory ranking, NOT true biological significance — with
# replicate data, swap this block out for a proper t-test/limma/DESeq2 call.
#
# IMPORTANT: with counts in the hundreds of millions, pnorm() underflows
# to exactly 0 for most points (-log10(0) = Inf). Compute in log-space
# with log.p = TRUE instead, which stays finite and preserves ranking.
df <- df %>%
  rowwise() %>%
  mutate(
    p_pool    = (EPI300_pCC1_ATF1_noAra + EPI300_pCC1_ATF1_iAAl) / (total_noAra + total_iAAl),
    se        = sqrt(p_pool * (1 - p_pool) * (1/total_noAra + 1/total_iAAl)),
    z         = (prop_iAAl - prop_noAra) / se,
    log_p_one_tail = pnorm(-abs(z), log.p = TRUE),           # natural log, one-tailed
    neglog10p_raw  = -(log_p_one_tail + log(2)) / log(10)    # two-tailed, base-10, finite
  ) %>%
  ungroup()

# ---- 4. Y-axis scale --------------------------------------------------
# -log10(p) reaches into the millions (a huge-N artifact, not real
# biological signal). Nothing is capped, clipped, or filtered: the axis is
# a pseudo-log scale (log-like for large values, linear near zero) whose
# breaks and limits are computed from the data, so the largest value
# always fits inside the panel.
df <- df %>% mutate(neglog10p = neglog10p_raw)

y_max <- max(df$neglog10p, na.rm = TRUE)

# Guard: if every point has neglog10p == 0 (e.g. Ara and noAra identical for
# all Names), log10(y_max) is -Inf and 0:ceiling(-Inf) would try to build an
# infinite vector. Fall back to a single break at 0 in that case.
if (!is.finite(y_max) || y_max <= 0) {
  y_breaks <- c(0)
} else {
  y_breaks <- c(0, 10^(0:ceiling(log10(y_max))))   # 0, 1, 10, ... past the max
}
y_limits <- c(0, max(y_breaks))

# X-axis view fixed to -8..8. Points beyond that (mostly Names found in only
# one condition) fall outside the visible window; they are not deleted or
# moved, the plot is just zoomed (coord_cartesian), so no warnings appear.
x_lim <- 6

# ---- 5. Flag significant hits ---------------------------------------
fc_cutoff <- 2       # |log2FC| threshold (dotted vertical lines at -2, 2)
p_cutoff  <- 1000     # -log10(p) threshold (dotted horizontal line at 1e3)

df <- df %>%
  mutate(significant = abs(log2FC) > fc_cutoff & neglog10p_raw > p_cutoff)

# Names that should always be labeled, regardless of significance
forced_labels <- c("1-Butanol, 3-methyl-")

# Only label a handful of the most extreme hits (by fold change, among
# significant points) so labels stay legible instead of overlapping.
top_labels <- df %>%
  filter(significant, !(Name %in% forced_labels), abs(log2FC) <= x_lim) %>%   # only label points inside the visible window
  distinct(Name, .keep_all = TRUE) %>%
  slice_max(abs(log2FC), n = 8, with_ties = FALSE) %>%
  pull(Name)

top_labels <- c(forced_labels, top_labels)

df <- df %>% mutate(label = ifelse(Name %in% top_labels, Name, NA))

# ---- 6. Plot ----------------------------------------------------------
p <- ggplot(df, aes(x = log2FC, y = neglog10p)) +
  geom_point(aes(color = significant), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey70")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
  geom_hline(yintercept = p_cutoff, linetype = "dotted") +
  geom_text_repel(
    data = ~ filter(.x, !is.na(label)),   # only the labelled points, so no NA rows are passed
    aes(label = label),
    size = 3.5,
    max.overlaps = Inf,
    segment.color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    force = 2,
    seed = 42
  ) +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(base = 10),
    breaks = y_breaks,
    limits = y_limits,
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_continuous(breaks = seq(-x_lim, x_lim, by = 2)) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  labs(
    title = "EPI300 pCC1-ATF1 Volatilome Response\n(no_arabinose vs iAAl)",
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~p~value~"(pseudo-log scale)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14)
  )

print(p)

ggsave("EPI300_pCC1_ATF1_noAraiAAL.png", p, width = 10, height = 8, dpi = 300)
