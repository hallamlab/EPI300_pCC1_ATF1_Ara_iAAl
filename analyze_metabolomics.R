## ============================================================================
## EPI300_pCC1_ATF1 metabolomics comparisons
## Heatmaps, top-abundant bar plots, and top up/down (log2FC) bar plots
## for: noAra vs Ara, noAra vs iAAl, noAra vs Ara+iAAl, Ara vs Ara+iAAl
## ============================================================================

## ---- packages --------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "pheatmap")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(pheatmap)

## ---- paths -------------------------------------------------------------
## Edit these two paths for your machine, or run this script from the folder
## that contains the CSVs, with "plots" as a subfolder for output.
DATA_DIR <- "data"
OUT_DIR  <- "plots"
dir.create(OUT_DIR, showWarnings = FALSE)

TOP_N_HEATMAP <- 30
TOP_N_ABUNDANT <- 20
TOP_N_UPDOWN <- 15

## ---- comparison definitions ------------------------------------------------
## Each entry: input csv, the two abundance columns to compare, and display labels
comparisons <- list(
  noAra_vs_Ara = list(
    file   = file.path(DATA_DIR, "EPI300_pCC1_ATF1_noAraAraOD.csv"),
    col1   = "EPI300_pCC1_ATF1_noAra",
    col2   = "EPI300_pCC1_ATF1_Ara",
    label1 = "noAra",
    label2 = "Ara"
  ),
  noAra_vs_iAAl = list(
    file   = file.path(DATA_DIR, "EPI300_pCC1_ATF1_noAraiAAlOD.csv"),
    col1   = "EPI300_pCC1_ATF1_noAra",
    col2   = "EPI300_pCC1_ATF1_iAAl",
    label1 = "noAra",
    label2 = "iAAl"
  ),
  noAra_vs_AraiAAl = list(
    file   = file.path(DATA_DIR, "EPI300_pCC1_ATF1_noAraAraiAAlOD.csv"),
    col1   = "EPI300_pCC1_ATF1_noAra",
    col2   = "EPI300_pCC1_ATF1_Ara_iAAl",
    label1 = "noAra",
    label2 = "Ara+iAAl"
  ),
  Ara_vs_AraiAAl = list(
    file   = file.path(DATA_DIR, "EPI300_pCC1_ATF1_AraAraiAAlOD.csv"),
    col1   = "EPI300_pCC1_ATF1_Ara",
    col2   = "EPI300_pCC1_ATF1_AraAraiAAl",
    label1 = "Ara",
    label2 = "Ara+iAAl"
  )
)

## ---- helpers ----------------------------------------------------------------

## Read + clean one comparison's CSV: drop empty trailing columns, coerce
## numerics, and collapse duplicate compound names by summing (co-eluting hits).
load_clean <- function(entry) {
  df <- read_csv(entry$file, show_col_types = FALSE)
  df <- df[, colSums(!is.na(df)) > 0]           # drop all-NA columns
  names(df)[1] <- "Name"
  df$Name <- str_trim(as.character(df$Name))
  df <- df[!is.na(df$Name) & df$Name != "", ]
  df[[entry$col1]] <- suppressWarnings(as.numeric(df[[entry$col1]]))
  df[[entry$col2]] <- suppressWarnings(as.numeric(df[[entry$col2]]))
  ## sum-with-min-count-1: NA only if every value in the group is NA
  sum_min1 <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  df <- df %>%
    group_by(Name) %>%
    summarise(
      !!entry$col1 := sum_min1(.data[[entry$col1]]),
      !!entry$col2 := sum_min1(.data[[entry$col2]]),
      .groups = "drop"
    ) %>%
    filter(!is.na(.data[[entry$col1]]) | !is.na(.data[[entry$col2]]))
  df
}

## Add log2 fold change (col2 / col1) using a shared pseudocount = half the
## smallest nonzero value observed across both columns (NA treated as 0).
add_fc <- function(df, entry) {
  c1 <- entry$col1; c2 <- entry$col2
  vals <- c(df[[c1]], df[[c2]])
  nonzero <- vals[!is.na(vals) & vals > 0]
  pseudo <- if (length(nonzero) > 0) min(nonzero) / 2 else 1
  a <- ifelse(is.na(df[[c1]]), 0, df[[c1]]) + pseudo
  b <- ifelse(is.na(df[[c2]]), 0, df[[c2]]) + pseudo
  df$log2FC <- log2(b / a)
  df$mean_abund <- (ifelse(is.na(df[[c1]]), 0, df[[c1]]) +
                       ifelse(is.na(df[[c2]]), 0, df[[c2]])) / 2
  attr(df, "pseudo") <- pseudo
  df
}

## Heatmap of the N most abundant compounds: log10(abundance + 1), two columns.
plot_heatmap <- function(df, entry, name) {
  c1 <- entry$col1; c2 <- entry$col2
  top <- df %>% arrange(desc(mean_abund)) %>% head(TOP_N_HEATMAP)
  mat <- log10(as.matrix(top[, c(c1, c2)]) %>% replace(is.na(.), 0) + 1)
  rownames(mat) <- str_trunc(top$Name, 45)
  colnames(mat) <- c(entry$label1, entry$label2)

  png(file.path(OUT_DIR, paste0(name, "_heatmap.png")),
      width = 7, height = max(6, 0.3 * nrow(mat)), units = "in", res = 200)
  pheatmap(mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("#440154", "#21908C", "#FDE725"))(100),
           main = paste0("Top ", nrow(mat), " most abundant compounds\n",
                          entry$label1, " vs ", entry$label2),
           fontsize_row = 8, angle_col = 0)
  dev.off()
}

## Top-N most abundant compounds (by mean of the two conditions), grouped bars.
plot_top_abundant <- function(df, entry, name) {
  c1 <- entry$col1; c2 <- entry$col2
  top <- df %>% arrange(desc(mean_abund)) %>% head(TOP_N_ABUNDANT)
  top$Name <- factor(top$Name, levels = rev(top$Name))

  long <- top %>%
    select(Name, all_of(c(c1, c2))) %>%
    pivot_longer(cols = all_of(c(c1, c2)), names_to = "Condition", values_to = "Abundance") %>%
    mutate(
      Condition = ifelse(Condition == c1, entry$label1, entry$label2),
      Abundance = log10(ifelse(is.na(Abundance), 0, Abundance) + 1)
    )

  p <- ggplot(long, aes(x = Name, y = Abundance, fill = Condition)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
    labs(title = paste0("Top ", nrow(top), " most abundant compounds (ranked by mean) — ",
                         entry$label1, " vs ", entry$label2),
         x = NULL, y = "log10(peak abundance + 1)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = c(0.85, 0.15))

  ggsave(file.path(OUT_DIR, paste0(name, "_top_abundant.png")), p,
         width = 9, height = max(6, 0.35 * nrow(top)), dpi = 200)
}

## Top-N up and top-N down compounds by log2 fold change, diverging bar plot.
plot_top_updown <- function(df, entry, name) {
  pseudo <- attr(df, "pseudo")
  d <- df %>% filter(!is.na(log2FC))
  up <- d %>% arrange(desc(log2FC)) %>% head(TOP_N_UPDOWN)
  down <- d %>% arrange(log2FC) %>% head(TOP_N_UPDOWN)
  combo <- bind_rows(up, down) %>%
    distinct(Name, .keep_all = TRUE) %>%
    arrange(log2FC) %>%
    mutate(
      label = str_trunc(Name, 45),
      label = factor(label, levels = label),
      direction = ifelse(log2FC > 0, "Up", "Down")
    )

  p <- ggplot(combo, aes(x = label, y = log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c(Up = "#d9534f", Down = "#4b7bec")) +
    labs(title = paste0("Top up/down compounds: ", entry$label2, " vs ", entry$label1),
         subtitle = sprintf("log2 fold change, pseudocount = %.1f", pseudo),
         x = NULL, y = paste0("log2FC (", entry$label2, " / ", entry$label1, ")")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")

  ggsave(file.path(OUT_DIR, paste0(name, "_top_updown.png")), p,
         width = 9, height = max(6, 0.3 * nrow(combo)), dpi = 200)
}

## ---- run all comparisons ----------------------------------------------------
summary_rows <- list()
for (nm in names(comparisons)) {
  entry <- comparisons[[nm]]
  df <- load_clean(entry)
  df <- add_fc(df, entry)

  plot_heatmap(df, entry, nm)
  plot_top_abundant(df, entry, nm)
  plot_top_updown(df, entry, nm)

  summary_rows[[nm]] <- data.frame(
    comparison = nm, n_compounds = nrow(df), pseudocount = attr(df, "pseudo")
  )
  message(sprintf("%s: %d compounds, pseudocount = %.3f", nm, nrow(df), attr(df, "pseudo")))
}

summary_df <- bind_rows(summary_rows)
print(summary_df)
write_csv(summary_df, file.path(OUT_DIR, "comparison_summary.csv"))
