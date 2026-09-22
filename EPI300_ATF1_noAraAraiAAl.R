# ============================================================
# Volcano plot: EPI300_pCC1_ATF1 (Ara_iAAl vs noAra)
# Styled to match Met_EPPCC1_EPFOS_iAA_Volcano.png
# ============================================================

library(tidyverse)   # dplyr, ggplot2, readr
library(ggrepel)     # non-overlapping labels with leader lines

# ---- 1. Load data -------------------------------------------------
# Expects a CSV with columns: Name, EPI300_pCC1_ATF1_noAra, EPI300_pCC1_ATF1_Ara_iAAl
df <- read_csv("EPI300_pCC1_ATF1_noAraAraiAAlOD.csv", show_col_types = FALSE)

# The export has trailing empty columns (from trailing commas in the CSV) and
# at least one row where an Excel formula error ("#VALUE!") landed in the
# Ara_iAAl column, which forces the whole column to be read as text (character)
# instead of numeric -- that's what made sum() fail with "invalid 'type'
# (character)". Coercing to numeric turns any such non-numeric cell into NA,
# and dropping the all-empty trailing columns cleans up the rest.
df <- df %>%
  select(Name, EPI300_pCC1_ATF1_noAra, EPI300_pCC1_ATF1_Ara_iAAl) %>%
  mutate(
    EPI300_pCC1_ATF1_noAra    = suppressWarnings(as.numeric(EPI300_pCC1_ATF1_noAra)),
    EPI300_pCC1_ATF1_Ara_iAAl = suppressWarnings(as.numeric(EPI300_pCC1_ATF1_Ara_iAAl))
  )

# ---- 2. Fold change -------------------------------------------------
# Blank cells and cells that failed numeric conversion (like "#VALUE!") are
# now NA. Treated as 0 counts here, with a small pseudocount added only
# inside the fold-change ratio, so no row needs to be dropped and log2(0)
# / division-by-zero can't occur. Every row in the file stays on the plot.
pseudo <- 0.5

df <- df %>%
  mutate(
    EPI300_pCC1_ATF1_noAra    = replace_na(EPI300_pCC1_ATF1_noAra, 0),
    EPI300_pCC1_ATF1_Ara_iAAl = replace_na(EPI300_pCC1_ATF1_Ara_iAAl, 0)
  )

total_noAra    <- sum(df$EPI300_pCC1_ATF1_noAra)
total_Ara_iAAl <- sum(df$EPI300_pCC1_ATF1_Ara_iAAl)

df <- df %>%
  mutate(
    prop_noAra    = EPI300_pCC1_ATF1_noAra    / total_noAra,
    prop_Ara_iAAl = EPI300_pCC1_ATF1_Ara_iAAl / total_Ara_iAAl,
    log2FC        = log2(((EPI300_pCC1_ATF1_Ara_iAAl + pseudo) / total_Ara_iAAl) /
                         ((EPI300_pCC1_ATF1_noAra    + pseudo) / total_noAra))
  )

# ---- 3. Significance ------------------------------------------------
# NOTE: this dataset has ONE measurement per Name per condition (no
# biological replicates), so a textbook p-value isn't really available.
# As a stand-in, this uses a two-proportion z-test of each Name's
# share of the total signal in Ara_iAAl vs noAra. Treat these p-values as a
# rough, exploratory ranking, NOT true biological significance — with
# replicate data, swap this block out for a proper t-test/limma/DESeq2 call.
#
# IMPORTANT: with counts in the hundreds of millions, pnorm() underflows
# to exactly 0 for most points (-log10(0) = Inf). Compute in log-space
# with log.p = TRUE instead, which stays finite and preserves ranking.
df <- df %>%
  rowwise() %>%
  mutate(
    p_pool    = (EPI300_pCC1_ATF1_noAra + EPI300_pCC1_ATF1_Ara_iAAl) / (total_noAra + total_Ara_iAAl),
    se        = sqrt(p_pool * (1 - p_pool) * (1/total_noAra + 1/total_Ara_iAAl)),
    # A row with 0 counts in both conditions gives se == 0 and a 0/0 = NaN
    # z-score, which would otherwise drop that point from the plot. Treat
    # "no signal in either condition" as z = 0 (no difference) instead.
    z         = ifelse(se == 0, 0, (prop_Ara_iAAl - prop_noAra) / se),
    log_p_one_tail = pnorm(-abs(z), log.p = TRUE),           # natural log, one-tailed
    neglog10p_raw  = -(log_p_one_tail + log(2)) / log(10)    # two-tailed, base-10, finite
  ) %>%
  ungroup()

# ---- 4. Y-axis scale --------------------------------------------------
# The raw values above span an enormous range (roughly 0 to several
# hundred thousand) — a huge-N artifact, not real biological signal
# (flagged earlier). A hard cap just pins everything to one flat line.
# Instead, plot on a "pseudo-log" scale: it behaves like log10 for large
# values but stays linear near zero, so points spread out across the
# panel instead of stacking at a ceiling. No capping/clipping of values
# needed — every point keeps its true relative position.
df <- df %>% mutate(neglog10p = neglog10p_raw)

y_max <- max(df$neglog10p, na.rm = TRUE)

# Guard: if every point ended up with neglog10p == 0, log10(y_max) is -Inf
# and 0:ceiling(-Inf) would try to build an infinite vector and crash.
if (!is.finite(y_max) || y_max <= 0) {
  y_breaks <- c(0)
} else {
  y_breaks <- c(0, 10^(0:ceiling(log10(y_max))))   # 0, 1, 10, ... past the max
}
y_limits <- c(0, max(y_breaks))

# X-axis view fixed to -8..8. Points beyond that fall outside the visible
# window; they are not deleted or moved, the plot is just zoomed
# (coord_cartesian below), so no "removed rows" warnings appear.
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
    title = "EPI300 pCC1-ATF1 Volatilome Response\n(no_arabinose vs arabinose and iAAl)",
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

ggsave("EPI300_pCC1_ATF1_noAraAraiAAl.png", p, width = 10, height = 8, dpi = 300)
