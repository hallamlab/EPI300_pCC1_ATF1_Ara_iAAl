# ============================================================
# Volcano plot: EPI300_pCC1_ATF1
# Ara_iAAL vs Ara
# ============================================================

library(tidyverse)
library(ggrepel)
library(scales)

# ---- 1. Load data ------------------------------------------------

df <- read_csv("EPI300_pCC1_ATF1_AraAraiAAlCla.csv",
               show_col_types = FALSE) %>%
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

top_labels <- df %>%
  filter(significant) %>%
  distinct(Name, .keep_all = TRUE) %>%
  slice_max(abs(log2FC),
            n = 15,
            with_ties = FALSE) %>%
  pull(Name)

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
    title = "EPI300 pCC1-ATF1 Volatilome Response\n(Ara vs Ara + iAAL)",
    x = expression(log[2] ~ fold ~ change),
    y = expression(-log[10] ~ p ~ value)
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16
    ),
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