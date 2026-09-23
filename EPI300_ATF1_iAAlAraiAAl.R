# ============================================================
# Volcano plot: EPI300_pCC1_ATF1
# Comparison: Ara + iAAl  vs  iAAl   (log2FC > 0 = higher with Ara + iAAl)
# ============================================================

library(tidyverse)
library(ggrepel)
library(scales)

# ---- 0. Settings (everything you might want to tweak) -----------

infile     <- "EPI300_pCC1_ATF1_iAAlAraiAAlOD.csv"
outfile    <- "EPI300_pCC1_ATF1_iAAlAraiAAL.png"

col_ref    <- "EPI300_pCC1_ATF1_iAAl"      # reference / denominator  (iAAl)
col_treat  <- "EPI300_pCC1_ATF1_AraiAAl"   # treatment / numerator    (Ara + iAAl)

pseudo     <- 0.5     # pseudocount, only used inside the fold-change ratio
fc_cutoff  <- 2       # |log2FC| threshold
p_cutoff   <- 1000    # -log10(p) threshold
x_lim      <- 6       # x-axis window; points beyond it are pinned to the edge
n_per_side <- 9       # top N labels for under-expressed (over-expressed: all are labelled)
max_chars  <- 34      # long compound names are shortened to this many characters

forced_labels <- c(
  "1-Butanol, 3-methyl-, acetate",
  "1-Butanol, 3-methyl-"
)

# This compound's label is placed separately (nudged left of the edge triangle);
# all labels get a line to their point.
line_labels <- "1-Butanol, 3-methyl-, acetate"

col_up   <- "#C0392B"   # higher with Ara + iAAl
col_down <- "#2471A3"   # lower with Ara + iAAl
col_ns   <- "grey75"

# ---- 1. Load data ------------------------------------------------

df <- read_csv(infile, show_col_types = FALSE) %>%
  select(Name, ref = all_of(col_ref), treat = all_of(col_treat)) %>%   # drops the empty extra column
  mutate(Name = str_trim(Name))

# Rows with no Name (blank rows at the end of the export) can't be plotted.
# They're reported instead of silently dropped.
n_unnamed <- sum(is.na(df$Name))
if (n_unnamed > 0) {
  message(n_unnamed, " rows have no Name and are excluded from the plot.")
}
df <- df %>% filter(!is.na(Name))

# Force the count columns to numeric (a stray "#VALUE!" from Excel makes
# read_csv import the whole column as text). Anything non-numeric or blank
# is treated as 0; the pseudocount below keeps log2(0) from happening.
df <- df %>%
  mutate(
    ref   = replace_na(suppressWarnings(as.numeric(ref)),   0),
    treat = replace_na(suppressWarnings(as.numeric(treat)), 0)
  )

# ---- 2. Fold change (on library-size-normalised proportions) -----

total_ref   <- sum(df$ref)
total_treat <- sum(df$treat)

df <- df %>%
  mutate(
    prop_ref   = ref   / total_ref,
    prop_treat = treat / total_treat,
    log2FC     = log2(((treat + pseudo) / total_treat) /
                      ((ref   + pseudo) / total_ref))
  )

# ---- 3. Significance (two-proportion z-test) ---------------------

df <- df %>%
  mutate(
    p_pool = (ref + treat) / (total_ref + total_treat),
    se     = sqrt(p_pool * (1 - p_pool) * (1 / total_ref + 1 / total_treat)),
    # 0 counts in both conditions -> se = 0 -> 0/0. Treat as "no difference".
    z      = if_else(se == 0, 0, (prop_treat - prop_ref) / se),
    # two-sided p on the log scale so huge |z| doesn't underflow to 0
    neglog10p = -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
  )

# ---- 4. Classification -------------------------------------------

df <- df %>%
  mutate(
    significant = abs(log2FC) > fc_cutoff & neglog10p > p_cutoff,
    group = case_when(
      significant & log2FC > 0 ~ "Higher with Ara + iAAl",
      significant & log2FC < 0 ~ "Lower with Ara + iAAl",
      TRUE                     ~ "Not significant"
    ),
    group = factor(group, levels = c("Higher with Ara + iAAl",
                                     "Lower with Ara + iAAl",
                                     "Not significant")),
    # Points beyond the x window are pinned to the edge and drawn as
    # triangles, so they stay visible instead of vanishing off the plot.
    off_scale = abs(log2FC) > x_lim,
    x_plot    = pmax(pmin(log2FC, x_lim), -x_lim)
  )

# Y-axis breaks (guarded for the all-zero / non-finite case)
y_max <- max(df$neglog10p, na.rm = TRUE)
if (!is.finite(y_max) || y_max <= 0) {
  y_breaks <- c(0)
} else {
  y_breaks <- c(0, 10^(0:floor(log10(y_max))))
}
y_limits <- c(0, if (is.finite(y_max) && y_max > 0) y_max * 1.5 else 1)

# ---- 5. Labels ---------------------------------------------------

# Forced labels always get a label; the rest are the most extreme
# significant hits on each side (by |log2FC|), excluding off-scale ones
# so labels don't pile up on the plot edge.
sig_df <- df %>% filter(significant, !(Name %in% forced_labels)) %>%
  distinct(Name, .keep_all = TRUE)

n_forced_under <- sum(df$Name %in% forced_labels & df$log2FC < 0)

# ALL significant over-expressed compounds are labelled (no cap).
top_over <- sig_df %>%
  filter(log2FC > 0) %>%
  pull(Name)

# Under-expressed: top N by |log2FC|, excluding off-scale ones so labels
# don't pile up on the plot edge.
top_under <- sig_df %>%
  filter(log2FC < 0, !off_scale) %>%
  slice_min(log2FC, n = max(n_per_side - n_forced_under, 0), with_ties = FALSE) %>%
  pull(Name)

missing_forced <- setdiff(forced_labels, df$Name)
if (length(missing_forced) > 0) {
  message("Forced label(s) not found in the data (check spelling): ",
          paste(missing_forced, collapse = "; "))
}

label_names <- c(forced_labels, top_over, top_under)

df <- df %>%
  mutate(label = if_else(Name %in% label_names,
                         str_trunc(Name, max_chars), NA_character_))

n_off <- sum(df$off_scale)

# ---- 6. Volcano plot ---------------------------------------------

# Compact axis labels (0, 1, 10, 100, 1K, 10K, 1M ...). Written by hand because
# label_number(scale_cut = ...) errors on the 0 break in some scales versions.
short_lab <- function(x) {
  x <- round(x)   # pseudo-log back-transform gives e.g. 999999.9999
  case_when(
    x >= 1e6 ~ paste0(x / 1e6, "M"),
    x >= 1e3 ~ paste0(x / 1e3, "K"),
    TRUE     ~ as.character(x)
  )
}

p <- ggplot(df, aes(x = x_plot, y = neglog10p)) +

  # threshold guides (drawn first so they sit behind the points)
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = p_cutoff,
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +

  # non-significant points underneath, significant on top
  geom_point(data = ~ filter(.x, group == "Not significant"),
             aes(colour = group, shape = off_scale),
             size = 2, alpha = 0.6) +
  geom_point(data = ~ filter(.x, group != "Not significant"),
             aes(colour = group, shape = off_scale),
             size = 2.4, alpha = 0.85) +

  # labels with a connector line to their point, no frame
  geom_text_repel(
    data = ~ filter(.x, !is.na(label), !(Name %in% line_labels)),
    aes(label = label),
    colour = "grey15",
    size = 3.2,
    max.overlaps = Inf,
    segment.colour = "grey35",    # connector "stick" from label to point
    segment.size = 0.3,
    min.segment.length = 0,       # always draw the stick, even for short ones
    box.padding = 0.5,
    point.padding = 0.3,
    force = 4,
    max.iter = 20000,
    seed = 42,
    show.legend = FALSE
  ) +

  # the acetate label (same style, nudged left so it clears the edge triangle)
  geom_text_repel(
    data = ~ filter(.x, !is.na(label), Name %in% line_labels),
    aes(label = label),
    colour = "grey15",
    size = 3.2,
    max.overlaps = Inf,
    segment.colour = "grey30",
    segment.size = 0.4,
    min.segment.length = 0,
    box.padding = 0.8,
    point.padding = 0.4,
    nudge_x = -1.5,
    seed = 42,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = c("Higher with Ara + iAAl" = col_up,
               "Lower with Ara + iAAl"  = col_down,
               "Not significant"        = col_ns),
    drop = FALSE
  ) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +

  scale_y_continuous(
    trans  = pseudo_log_trans(base = 10),
    breaks = y_breaks,
    labels = short_lab,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_continuous(breaks = seq(-x_lim, x_lim, by = 2)) +
  # coord_cartesian only zooms -- unlike scale limits it never turns points or
  # labels into NA, so nothing (e.g. the acetate label) can be dropped.
  coord_cartesian(xlim = c(-x_lim - 0.3, x_lim + 0.3),
                  ylim = y_limits) +

  labs(
    title    = "EPI300 pCC1-ATF1 Volatilome Response",
    subtitle = "Ara + iAAl vs iAAl",
    x        = expression(log[2] ~ fold ~ change),
    y        = expression(-log[10] ~ p ~ value ~ "(pseudo-log scale)"),
    colour   = NULL,
    caption  = if (n_off > 0)
      paste0("Triangles: ", n_off, " compound(s) with |log2FC| > ", x_lim,
             " pinned to the plot edge.") else NULL
  ) +

  theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle    = element_text(hjust = 0.5),
    plot.caption     = element_text(colour = "grey40", size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    legend.position  = "none"
  )

print(p)

ggsave(outfile, p, width = 11, height = 8.5, dpi = 300)

# Quick summary in the console
print(table(df$group))
