# ============================================================
# Volcano plot: EPI300_pCC1_ATF1
# Ara_iAAL vs Ara
# ============================================================

library(tidyverse)
library(ggrepel)
library(scales)

# ---- 1. Load data ------------------------------------------------

df <- read_csv("EPI300_pCC1_ATF1_AraAraiAAlOD.csv",
               show_col_types = FALSE)

# The export has a handful of fully-blank trailing rows (no Name and no
# counts at all) -- those aren't real data points, just artifacts of the
# export, so they're dropped here (nothing with an actual Name is removed).
df <- df %>% filter(!is.na(Name))

# Force the count columns to numeric. A stray non-numeric entry
# (e.g. "#VALUE!", an Excel formula error, instead of a peak area)
# will make read_csv() import the WHOLE column as character, which
# then breaks sum()/arithmetic below. suppressWarnings() here is
# intentional: any cell that can't convert becomes NA.
df <- df %>%
  mutate(
    EPI300_pCC1_ATF1_Ara      = suppressWarnings(as.numeric(EPI300_pCC1_ATF1_Ara)),
    EPI300_pCC1_ATF1_AraAraiAAl = suppressWarnings(as.numeric(EPI300_pCC1_ATF1_AraAraiAAl))
  )

# Blank cells and cells that failed numeric conversion (like "#VALUE!") are
# now NA. Rather than dropping those Names (which would remove them from
# the plot), they're treated as 0 counts -- a pseudocount is added only
# inside the fold-change ratio below so log2(0) / division-by-zero can't
# happen. Every real Name stays on the plot.
df <- df %>%
  mutate(
    EPI300_pCC1_ATF1_Ara        = replace_na(EPI300_pCC1_ATF1_Ara, 0),
    EPI300_pCC1_ATF1_AraAraiAAl = replace_na(EPI300_pCC1_ATF1_AraAraiAAl, 0)
  )

# ---- 2. Fold Change ---------------------------------------------

pseudo <- 0.5

total_Ara  <- sum(df$EPI300_pCC1_ATF1_Ara)
total_iAAL <- sum(df$EPI300_pCC1_ATF1_AraAraiAAl)

df <- df %>%
  mutate(
    prop_Ara     = EPI300_pCC1_ATF1_Ara / total_Ara,
    prop_AraiAAL = EPI300_pCC1_ATF1_AraAraiAAl / total_iAAL,
    log2FC       = log2(((EPI300_pCC1_ATF1_AraAraiAAl + pseudo) / total_iAAL) /
                        ((EPI300_pCC1_ATF1_Ara         + pseudo) / total_Ara))
  )

# ---- 3. Significance --------------------------------------------

df <- df %>%
  rowwise() %>%
  mutate(
    p_pool =
      (EPI300_pCC1_ATF1_Ara +
         EPI300_pCC1_ATF1_AraAraiAAl) /
      (total_Ara + total_iAAL),
    
    se =
      sqrt(
        p_pool * (1 - p_pool) *
          (1 / total_Ara + 1 / total_iAAL)
      ),
    
    # A row with 0 counts in both conditions gives se == 0 and a 0/0 = NaN
    # z-score, which would otherwise drop that point from the plot. Treat
    # "no signal in either condition" as z = 0 (no difference) instead.
    z = ifelse(se == 0, 0, (prop_AraiAAL - prop_Ara) / se),
    
    log_p_one_tail =
      pnorm(-abs(z), log.p = TRUE),
    
    neglog10p =
      -(log_p_one_tail + log(2)) / log(10)
  ) %>%
  ungroup()

# ---- 4. Significance thresholds ---------------------------------

fc_cutoff <- 2
p_cutoff  <- 1000

df <- df %>%
  mutate(
    significant =
      abs(log2FC) > fc_cutoff &
      neglog10p > p_cutoff
  )

# Y-axis breaks computed from the data (guarded against the all-zero /
# non-finite case, where log10(y_max) would be -Inf and 0:ceiling(-Inf)
# would try to build an infinite vector).
y_max <- max(df$neglog10p, na.rm = TRUE)
if (!is.finite(y_max) || y_max <= 0) {
  y_breaks <- c(0)
} else {
  y_breaks <- c(0, 10^(0:ceiling(log10(y_max))))
}
y_limits <- c(0, max(y_breaks))

# X-axis view fixed to -6..6, matching your other volcano plots. This is a
# zoom (coord_cartesian, applied in the plot below), not a filter, so
# points beyond the window stay in the data and are never dropped.
x_lim <- 6

# ---- 5. Labels --------------------------------------------------

# Compounds that should always be labeled, regardless of significance
forced_labels <- c(
  "1-Butanol, 3-methyl-, acetate",
  "1-Butanol, 3-methyl-"
)

n_per_side <- 9  # top N over-expressed + top N under-expressed

# Figure out which side each forced label belongs to, from its own
# log2FC sign (not just among the "significant" points), so it still
# gets labeled even if it falls just short of the significance cutoff
forced_sign <- df %>%
  filter(Name %in% forced_labels) %>%
  distinct(Name, .keep_all = TRUE) %>%
  select(Name, log2FC)

forced_over  <- forced_sign %>% filter(log2FC > 0) %>% pull(Name)
forced_under <- forced_sign %>% filter(log2FC < 0) %>% pull(Name)

sig_df <- df %>%
  filter(significant) %>%
  distinct(Name, .keep_all = TRUE)

# Over-expressed (log2FC > 0): reserve slots for any forced labels on
# this side, then fill the rest with the largest positive log2FC values
top_over <- sig_df %>%
  filter(log2FC > 0, !(Name %in% forced_labels), log2FC <= x_lim) %>%
  slice_max(log2FC,
            n = max(n_per_side - length(forced_over), 0),
            with_ties = FALSE) %>%
  pull(Name)

# Under-expressed (log2FC < 0): same idea, using the most negative values
top_under <- sig_df %>%
  filter(log2FC < 0, !(Name %in% forced_labels), log2FC >= -x_lim) %>%
  slice_min(log2FC,
            n = max(n_per_side - length(forced_under), 0),
            with_ties = FALSE) %>%
  pull(Name)

top_labels <- c(forced_over, top_over, forced_under, top_under)

df <- df %>%
  mutate(
    label = ifelse(Name %in% top_labels,
                   Name,
                   NA)
  )

# ---- 6. Volcano Plot --------------------------------------------

p <- ggplot(df,
            aes(x = log2FC,
                y = neglog10p)) +
  
  geom_point(
    aes(color = significant),
    size = 2
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "red",
      "FALSE" = "grey70"
    )
  ) +
  
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    linetype = "dotted"
  ) +
  
  geom_hline(
    yintercept = p_cutoff,
    linetype = "dotted"
  ) +
  
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
    trans = pseudo_log_trans(base = 10),
    breaks = y_breaks,
    limits = y_limits,
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  scale_x_continuous(breaks = seq(-x_lim, x_lim, by = 2)) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  
  labs(
    title = "EPI300 pCC1-ATF1 Volatilome Response",
    subtitle = "(arabinose vs arabinose and iAAl)",
    x = expression(log[2] ~ fold ~ change),
    y = expression(-log[10] ~ p ~ value ~ "(pseudo-log scale)")
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16
    ),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p)

ggsave("EPI300_ATF1_AraAraiAAL.png",  p,width = 10,  height = 8,  dpi = 300)
