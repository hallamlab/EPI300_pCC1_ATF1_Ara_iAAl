# ============================================================
# Volcano plot: EPI300_pCC1_ATF1
# Ara_iAAL vs Ara
# ============================================================

library(tidyverse)
library(ggrepel)
library(scales)

# ---- 1. Load data ------------------------------------------------

df <- read_csv("EPI300_pCC1_ATF1_AraAraiAAlCla.csv",
               show_col_types = FALSE)

# Force the count columns to numeric. A stray non-numeric entry
# (e.g. "Find by Chromatogram Deconvolution" instead of a peak area)
# will make read_csv() import the WHOLE column as character, which
# then breaks sum()/arithmetic below. suppressWarnings() here is
# intentional: any cell that can't convert becomes NA, and gets
# reported and dropped in the next step.
df <- df %>%
  mutate(
    EPI300_pCC1_ATF1_Ara      = suppressWarnings(as.numeric(EPI300_pCC1_ATF1_Ara)),
    EPI300_pCC1_ATF1_Ara_iAAL = suppressWarnings(as.numeric(EPI300_pCC1_ATF1_Ara_iAAL))
  )

dropped <- df %>%
  filter(is.na(EPI300_pCC1_ATF1_Ara) | is.na(EPI300_pCC1_ATF1_Ara_iAAL))

if (nrow(dropped) > 0) {
  message(sprintf("Dropping %d row(s) with missing/non-numeric values:", nrow(dropped)))
  print(dropped %>% select(Name, EPI300_pCC1_ATF1_Ara, EPI300_pCC1_ATF1_Ara_iAAL),
        n = nrow(dropped))
}

df <- df %>%
  filter(
    !is.na(EPI300_pCC1_ATF1_Ara),
    !is.na(EPI300_pCC1_ATF1_Ara_iAAL)
  )

# ---- 2. Fold Change ---------------------------------------------

total_Ara <- sum(df$EPI300_pCC1_ATF1_Ara)
total_iAAL <- sum(df$EPI300_pCC1_ATF1_Ara_iAAL)

df <- df %>%
  mutate(
    prop_Ara  = EPI300_pCC1_ATF1_Ara / total_Ara,
    prop_iAAL = EPI300_pCC1_ATF1_Ara_iAAL / total_iAAL,
    log2FC    = log2(prop_iAAL / prop_Ara)
  )

# ---- 3. Significance --------------------------------------------

df <- df %>%
  rowwise() %>%
  mutate(
    p_pool =
      (EPI300_pCC1_ATF1_Ara +
         EPI300_pCC1_ATF1_Ara_iAAL) /
      (total_Ara + total_iAAL),
    
    se =
      sqrt(
        p_pool * (1 - p_pool) *
          (1 / total_Ara + 1 / total_iAAL)
      ),
    
    z = (prop_iAAL - prop_Ara) / se,
    
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

# ---- 5. Labels --------------------------------------------------

# Compounds that should always be labeled, regardless of significance
forced_labels <- c(
  "1-Butanol, 3-methyl-, acetate",
  "1-Butanol, 3-methyl-"
)

n_per_side <- 4  # top N over-expressed + top N under-expressed

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
  filter(log2FC > 0, !(Name %in% forced_labels)) %>%
  slice_max(log2FC,
            n = max(n_per_side - length(forced_over), 0),
            with_ties = FALSE) %>%
  pull(Name)

# Under-expressed (log2FC < 0): same idea, using the most negative values
top_under <- sig_df %>%
  filter(log2FC < 0, !(Name %in% forced_labels)) %>%
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
    breaks = c(
      0, 2, 10, 100,
      1000, 10000, 100000
    )
  ) +
  
  labs(
    title = "EPI300 pCC1-ATF1 Volatilome Response",
    subtitle = "(no_arabinose vs arabinose and iAAI)",
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

ggsave(
  "EPI300_pCC1_ATF1_Ara_vs_Ara_iAAL_Volcano.png",
  p,
  width = 10,
  height = 8,
  dpi = 300
)