# ============================================================
# Volcano plot: EPI300_pCC1_ATF1  (no_arabinose vs iAAL)
# ============================================================

library(tidyverse)
library(ggrepel)
library(scales)

# ---- 1. Load data ------------------------------------------------
df <- read_csv("EPI300_pCC1_ATF1_noAraiAAlODCla.csv",
               show_col_types = FALSE) %>%
  select(Name, noAra = EPI300_pCC1_ATF1_noAra, iAAL = EPI300_pCC1_ATF1_iAAL) %>%
  mutate(across(c(noAra, iAAL), ~ suppressWarnings(as.numeric(.x))))

n_start <- nrow(df)

df <- df %>%
  filter(!is.na(Name), Name != "0",          # empty / placeholder names
         !is.na(noAra), !is.na(iAAL),        # non-numeric entries
         noAra > 0, iAAL > 0)                # keep only compounds detected in BOTH samples

message(sprintf("Kept %d of %d rows (removed empty, placeholder '0', and compounds detected in only one sample).",
                nrow(df), n_start))

# ---- 2. Fold change ---------------------------------------------
total_noAra <- sum(df$noAra)
total_iAAL  <- sum(df$iAAL)

df <- df %>%
  mutate(
    prop_noAra = noAra / total_noAra,
    prop_iAAL  = iAAL  / total_iAAL,
    log2FC     = log2(prop_iAAL / prop_noAra)
  )

# ---- 3. Significance --------------------------------------------
# One measurement per compound per sample (no replicates), so this is
# an exploratory two-proportion z-test, NOT a true biological p-value.
# Computed in log-space to avoid underflow to 0.
df <- df %>%
  mutate(
    p_pool    = (noAra + iAAL) / (total_noAra + total_iAAL),
    se        = sqrt(p_pool * (1 - p_pool) * (1 / total_noAra + 1 / total_iAAL)),
    z         = (prop_iAAL - prop_noAra) / se,
    neglog10p = -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
  )

fc_cutoff <- 2
p_cutoff  <- 1000

df <- df %>%
  mutate(significant = abs(log2FC) > fc_cutoff & neglog10p > p_cutoff,
         direction   = case_when(significant & log2FC > 0 ~ "Higher in iAAL",
                                 significant & log2FC < 0 ~ "Higher in no_Ara",
                                 TRUE                     ~ "Not significant"))

# ---- 4. Labels --------------------------------------------------
forced_labels <- c("1-Butanol, 3-methyl-")   # always label these
n_per_side    <- 4                           # top N per direction

top_over  <- df %>% filter(significant, log2FC > 0) %>%
  slice_max(neglog10p, n = n_per_side, with_ties = FALSE) %>% pull(Name)
top_under <- df %>% filter(significant, log2FC < 0) %>%
  slice_max(neglog10p, n = n_per_side, with_ties = FALSE) %>% pull(Name)

df <- df %>%
  mutate(label = if_else(Name %in% c(forced_labels, top_over, top_under),
                         str_wrap(str_trunc(Name, 45), width = 22),
                         NA_character_))

# Labeled significant points get pushed outward, into the empty space
# beside the plot, so text doesn't sit on top of the grey cloud.
label_df <- df %>%
  filter(!is.na(label)) %>%
  mutate(nudge_x = case_when(
    abs(log2FC) <= fc_cutoff                 ~ 0,                    # forced label near centre
    abs(log2FC) > max(abs(log2FC)) - 1.5     ~ -sign(log2FC) * 2,    # near the edge: push inward
    TRUE                                     ~  sign(log2FC) * 2     # otherwise push outward
  ))

# ---- 5. Axes ----------------------------------------------------
x_lim <- 2 * ceiling((max(abs(df$log2FC)) + 1) / 2)   # symmetric, nothing clipped
y_brk <- c(0, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7)
y_brk <- y_brk[y_brk <= max(df$neglog10p) * 10]
y_lab <- parse(text = c("0", "10", "10^2", "10^3", "10^4",
                        "10^5", "10^6", "10^7"))[seq_along(y_brk)]

# ---- 6. Plot ----------------------------------------------------
p <- ggplot(df, aes(x = log2FC, y = neglog10p)) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
  geom_hline(yintercept = p_cutoff, linetype = "dotted") +
  geom_point(aes(color = direction),
             size = 2.4, alpha = 0.85) +
  scale_color_manual(
    values = c("Higher in iAAL"   = "red",
               "Higher in no_Ara" = "blue",
               "Not significant"  = "grey70"),
    breaks = c("Higher in iAAL", "Higher in no_Ara")
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = label),
    nudge_x = label_df$nudge_x,
    xlim = c(-x_lim, x_lim),       # keep labels inside the panel
    size = 3.2, lineheight = 0.9,
    max.overlaps = Inf,
    segment.color = "black", segment.size = 0.3,
    box.padding = 0.9, point.padding = 0.5,
    min.segment.length = 0, force = 4, seed = 42
  ) +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = seq(-x_lim, x_lim, by = 2),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                     breaks = y_brk, labels = y_lab,
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title    = "EPI300 pCC1-ATF1 Volatilome Response",
    subtitle = "(no_arabinose vs iAAL)",
    x        = expression(log[2]~fold~change~"(iAAL / no_Ara)"),
    y        = expression(-log[10]~p~value~"(pseudo-log scale)"),
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 16),
    plot.subtitle    = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

print(p)

ggsave("EPI300_ATF1_noAraiAAL.png", p, width = 10, height = 8, dpi = 300)
