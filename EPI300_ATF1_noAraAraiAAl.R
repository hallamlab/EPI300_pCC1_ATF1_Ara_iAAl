# ============================================================
# Volcano plot: EPI300_pCC1_ATF1 (Ara_iAAL vs noAra)
# Styled to match Met_EPPCC1_EPFOS_iAA_Volcano.png
# ============================================================

library(tidyverse)   # dplyr, ggplot2, readr
library(ggrepel)     # non-overlapping labels with leader lines

# ---- 1. Load data -------------------------------------------------
# Expects a CSV with columns: Compound, EPI300_pCC1_ATF1_noAra, EPI300_pCC1_ATF1_Ara_iAAL
df <- read_csv("EPI300_pCC1_ATF1_noAraAraiAAlCla.csv", show_col_types = FALSE) %>%
  filter(Compound != "0", !is.na(EPI300_pCC1_ATF1_noAra), !is.na(EPI300_pCC1_ATF1_Ara_iAAL))

# ---- 2. Fold change -------------------------------------------------
total_noAra <- sum(df$EPI300_pCC1_ATF1_noAra)
total_AraiAAL   <- sum(df$EPI300_pCC1_ATF1_Ara_iAAL)

df <- df %>%
  mutate(
    prop_noAra   = EPI300_pCC1_ATF1_noAra   / total_noAra,
    prop_AraiAAL = EPI300_pCC1_ATF1_Ara_iAAL / total_AraiAAL,
    log2FC       = log2(prop_AraiAAL / prop_noAra)
  )

# ---- 3. Significance ------------------------------------------------
# NOTE: this dataset has ONE measurement per compound per condition (no
# biological replicates), so a textbook p-value isn't really available.
# As a stand-in, this uses a two-proportion z-test of each compound's
# share of the total signal in Ara_iAAL vs noAra. Treat these p-values as a
# rough, exploratory ranking, NOT true biological significance — with
# replicate data, swap this block out for a proper t-test/limma/DESeq2 call.
#
# IMPORTANT: with counts in the hundreds of millions, pnorm() underflows
# to exactly 0 for most points (-log10(0) = Inf). Compute in log-space
# with log.p = TRUE instead, which stays finite and preserves ranking.
df <- df %>%
  rowwise() %>%
  mutate(
    p_pool    = (EPI300_pCC1_ATF1_noAra + EPI300_pCC1_ATF1_Ara_iAAL) / (total_noAra + total_AraiAAL),
    se        = sqrt(p_pool * (1 - p_pool) * (1/total_noAra + 1/total_AraiAAL)),
    z         = (prop_AraiAAL - prop_noAra) / se,
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

# ---- 5. Flag significant hits ---------------------------------------
fc_cutoff <- 2       # |log2FC| threshold (dotted vertical lines at -2, 2)
p_cutoff  <- 1000     # -log10(p) threshold (dotted horizontal line at 1e3)

df <- df %>%
  mutate(significant = abs(log2FC) > fc_cutoff & neglog10p_raw > p_cutoff)

# Only label a handful of the most extreme hits (by fold change, among
# significant points) so labels stay legible instead of overlapping.
top_labels <- df %>%
  filter(significant) %>%
  distinct(Compound, .keep_all = TRUE) %>%
  slice_max(abs(log2FC), n = 8, with_ties = FALSE) %>%
  pull(Compound)

df <- df %>% mutate(label = ifelse(Compound %in% top_labels, Compound, NA))

# ---- 6. Plot ----------------------------------------------------------
p <- ggplot(df, aes(x = log2FC, y = neglog10p)) +
  geom_point(aes(color = significant), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey70")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
  geom_hline(yintercept = p_cutoff, linetype = "dotted") +
  geom_text_repel(
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
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 2, 10, 100, 1000, 10000, 100000)
  ) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 2)) +
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

ggsave("EPI300_ATF1_noAraAraiAAl.png", p, width = 10, height = 8, dpi = 300)
