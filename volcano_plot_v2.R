# ============================================================
# Volcano-style plot: EPI300 pCC1-ATF1 Volatilome Response
# (no_arabinose vs arabinose and iAAI)
#
# Built for a CSV like: Name, EPI300_pCC1_ATF1_Ara, EPI300_pCC1_ATF1_Ara_iAAL
# i.e. ONE raw peak-area value per compound per condition -- no
# replicates, so no real statistical p-value can be computed.
#
# IMPORTANT: because there's no replicate data, the p-value on the
# y-axis is NOT from a real biological-replicate statistical test.
# It's computed with a Poisson two-sample test (stats::poisson.test),
# treating each peak area as a Poisson-distributed count -- the same
# logic behind the classic Audic-Claverie test for comparing single
# counts. It's a reasonable proxy for ranking/screening compounds,
# but it should NOT be presented as equivalent to a p-value from
# actual biological replicates. If you obtain replicate data later,
# swap in a real t-test (see the very first version of this script).
# ============================================================

library(ggplot2)
library(ggrepel)
library(scales)   # for pseudo_log_trans()

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

df <- read.csv("EPI300_pCC1_ATF1_AraAraiAAlCla.csv", stringsAsFactors = FALSE)

# Your actual column names -- change these if you rename the file
COMPOUND_COL <- "Name"
COND1_COL    <- "EPI300_pCC1_ATF1_Ara"       # e.g. arabinose
COND2_COL    <- "EPI300_pCC1_ATF1_Ara_iAAL"  # e.g. arabinose + iAAL

missing_cols <- setdiff(c(COMPOUND_COL, COND1_COL, COND2_COL), names(df))
if (length(missing_cols) > 0) {
  stop(
    "These columns were not found in your CSV: ",
    paste(missing_cols, collapse = ", "),
    "\nColumns actually present: ", paste(names(df), collapse = ", ")
  )
}

names(df)[names(df) == COMPOUND_COL] <- "Compound"
names(df)[names(df) == COND1_COL]    <- "cond1"
names(df)[names(df) == COND2_COL]    <- "cond2"

# ------------------------------------------------------------
# 2. Clean up: force numeric, drop rows with an unusable value
#    (e.g. a stray text entry like "Find by Chromatogram Deconvolution")
# ------------------------------------------------------------

raw_cond1 <- df$cond1
raw_cond2 <- df$cond2

df$cond1 <- suppressWarnings(as.numeric(df$cond1))
df$cond2 <- suppressWarnings(as.numeric(df$cond2))

bad_cond1 <- is.na(df$cond1) & !is.na(raw_cond1) & raw_cond1 != ""
bad_cond2 <- is.na(df$cond2) & !is.na(raw_cond2) & raw_cond2 != ""

if (any(bad_cond1) || any(bad_cond2)) {
  message("Dropping rows with non-numeric values, e.g.:")
  print(df[bad_cond1 | bad_cond2, c("Compound", "cond1", "cond2")],
        row.names = FALSE)
}

# Compounds only detected in ONE condition (blank in the other) get
# treated as missing rather than zero, since a blank here means "not
# found," not "found at zero abundance." Change this if you'd rather
# treat blanks as true zeros.
n_before <- nrow(df)
df <- df[!(is.na(df$cond1) & is.na(df$cond2)), ]
n_after <- nrow(df)
if (n_after < n_before) {
  message(sprintf("Dropped %d row(s) with no usable value in either condition.",
                   n_before - n_after))
}

# ------------------------------------------------------------
# 3. Compute log2 fold change and a p value
# ------------------------------------------------------------
# NOTE ON THE P VALUE: with only one measurement per compound per
# condition (no replicates), a true biological-replicate p value is
# not possible. As a stand-in, this uses a Poisson two-sample test
# (stats::poisson.test) -- the same idea behind the classic
# Audic-Claverie test used for comparing single-library sequencing
# counts. It treats each peak area as a Poisson-distributed count and
# asks how likely a difference this large is if the true rates were
# equal. This is a reasonable, commonly-used proxy, but it is NOT
# equivalent to a p value from real biological replicates -- treat
# these as a ranking/screening tool, not a rigorous significance
# test, and confirm any hits of interest with replicated experiments.

PSEUDOCOUNT <- 1  # avoids log2(0) and division by zero

c1 <- ifelse(is.na(df$cond1), 0, df$cond1) + PSEUDOCOUNT
c2 <- ifelse(is.na(df$cond2), 0, df$cond2) + PSEUDOCOUNT

df$log2FC <- log2(c2 / c1)

poisson_pvalue <- function(x, y) {
  tryCatch({
    stats::poisson.test(c(round(x), round(y)))$p.value
  }, error = function(e) NA_real_)
}

df$pvalue <- mapply(poisson_pvalue, c1, c2)
df$neglog10p <- -log10(df$pvalue)

# ------------------------------------------------------------
# 4. Thresholds and significance flag
# ------------------------------------------------------------

LOG2FC_THRESHOLD <- 2       # vertical dashed lines at +/- this value
PVAL_THRESHOLD   <- 1e-3    # horizontal dashed line at this p value

df$significant <- abs(df$log2FC) > LOG2FC_THRESHOLD &
                   df$pvalue < PVAL_THRESHOLD

# ------------------------------------------------------------
# 5. Choose which points to label
# ------------------------------------------------------------

n_label <- 7

label_df <- df[df$significant, ]
label_df <- label_df[order(-abs(label_df$log2FC)), ]
label_df <- head(label_df, n_label)

df$label <- ifelse(df$Compound %in% label_df$Compound, df$Compound, NA)

# ------------------------------------------------------------
# 6. Build the plot
# ------------------------------------------------------------

p <- ggplot(df, aes(x = log2FC, y = neglog10p)) +
  geom_point(aes(color = significant), size = 2.5, alpha = 0.9) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey60"),
                      guide = "none") +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
             linetype = "dotted", color = "black") +
  geom_hline(yintercept = -log10(PVAL_THRESHOLD),
             linetype = "dotted", color = "black") +
  geom_text_repel(aes(label = label),
                   size = 3.5,
                   max.overlaps = Inf,
                   box.padding = 0.5,
                   segment.color = "grey40",
                   na.rm = TRUE) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                      breaks = c(0, 2, 10, 1e2, 1e3, 1e4, 1e5),
                      labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "EPI300 pCC1-ATF1 Volatilome Response",
    subtitle = "(no_arabinose vs arabinose and iAAI)",
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~p~value~"(pseudo-log scale)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p)

# ------------------------------------------------------------
# 7. Save to file
# ------------------------------------------------------------

ggsave("volcano_plot.png", plot = p, width = 10, height = 8, dpi = 300)
