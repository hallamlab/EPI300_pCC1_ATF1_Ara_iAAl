# ============================================================
# Extra visualizations: EPI300_pCC1_ATF1 volatilome
# noAra, Ara, iAAL, Ara + iAAL  --  all 4 pairwise comparisons
#
# For EACH comparison this makes 5 PNG figures:
#   1. Scatter plot        (A vs B, log10 axes)
#   2. MA plot              (log2FC vs average intensity)
#   3. Ranked bar chart     (top compounds up / down)
#   4. Heatmap              (top compounds by |log2FC|)
#   5. Top abundant bars    (most abundant compounds per condition)
# Plus two combined overview figures (scatter + MA, all 4 comparisons
# side by side) so you can eyeball everything at once.
#
# Each comparison file has: column 1 = compound name, column 2 = culture A,
# column 3 = culture B (OD-normalized intensities; 0 / blank = not detected).
# Columns are read BY POSITION, so header spelling doesn't matter.
# ============================================================

library(tidyverse)
library(ggrepel)
library(scales)

# ---- 1. What to compare -------------------------------------------
comparisons <- tribble(
  ~comparison,             ~file,                                    ~A,      ~B,
  "noAra vs Ara",          "EPI300_pCC1_ATF1_noAraAraOD.csv",         "noAra", "Ara",
  "noAra vs iAAL",         "EPI300_pCC1_ATF1_noAraiAAlOD.csv",        "noAra", "iAAL",
  "noAra vs Ara + iAAL",   "EPI300_pCC1_ATF1_noAraAraiAAlOD.csv",     "noAra", "Ara + iAAL",
  "iAAL vs Ara + iAAL",    "EPI300_pCC1_ATF1_iAAlAraiAAlOD.csv",      "iAAL",  "Ara + iAAL"
)

# Skip any comparison whose CSV isn't in the working directory (with a message),
# so one missing file doesn't stop the others from being plotted.
missing_files <- comparisons %>% filter(!file.exists(file))
if (nrow(missing_files) > 0)
  message("Skipping (file not found in ", getwd(), "): ",
          paste(missing_files$file, collapse = ", "))
comparisons <- comparisons %>% filter(file.exists(file))
if (nrow(comparisons) == 0) stop("None of the input CSV files were found in ", getwd())

culture_cols <- c("noAra" = "blue", "Ara" = "red",
                  "iAAL" = "darkgreen", "Ara + iAAL" = "purple")

fc_cutoff <- 2   # log2FC threshold (= 4-fold) used to call "Higher in ..."

# tidy filename tag for a comparison label, e.g. "iAAl vs Ara + iAAL" -> "iAAl_vs_Ara_iAAL"
bad_names <- setdiff(c(comparisons$A, comparisons$B), names(culture_cols))
if (length(bad_names) > 0)
  stop("Culture name(s) not in culture_cols (check spelling/case): ",
       paste(bad_names, collapse = ", "))

tag_of <- function(x) x %>% str_replace_all(" \\+ ", "_") %>% str_replace_all(" vs ", "_vs_")

base_theme <- theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())

# ---- 2. Load & prepare one comparison -------------------------------
# A blank cell means "not detected". A cell that is non-blank but still
# fails to parse (e.g. an Excel "#VALUE!" error) is a real problem and
# gets dropped with a message, not silently zeroed. Only compounds
# detected (> 0) in BOTH cultures are kept, since fold change needs both.
read_pair <- function(comparison, file, A_name, B_name) {
  raw <- read_csv(file, col_types = cols(.default = col_character()),
                  show_col_types = FALSE, name_repair = "minimal")[, 1:3]
  names(raw) <- c("Name", "A_raw", "B_raw")

  d <- raw %>%
    filter(!is.na(Name), Name != "0") %>%
    mutate(A_num = suppressWarnings(as.numeric(A_raw)),
           B_num = suppressWarnings(as.numeric(B_raw)),
           is_bad = (is.na(A_num) & !is.na(A_raw) & A_raw != "") |
                    (is.na(B_num) & !is.na(B_raw) & B_raw != ""))

  n_bad <- sum(d$is_bad)
  if (n_bad > 0)
    message(sprintf("[%s] dropped %d row(s) with non-numeric (non-blank) values.",
                    comparison, n_bad))

  d %>%
    filter(!is_bad) %>%
    mutate(Name = make.unique(Name, sep = " #"),      # repeated names = separate peaks
           A = replace_na(A_num, 0),
           B = replace_na(B_num, 0)) %>%
    filter(A > 0, B > 0) %>%                          # detected in BOTH cultures
    mutate(comparison = comparison, A_name = A_name, B_name = B_name) %>%
    select(comparison, Name, A, B, A_name, B_name)
}

all_calc <- pmap_dfr(comparisons, function(comparison, file, A, B)
                       read_pair(comparison, file, A, B)) %>%
  mutate(comparison = factor(comparison, levels = comparisons$comparison)) %>%
  group_by(comparison) %>%
  mutate(log2FC    = log2((B / sum(B)) / (A / sum(A))),   # share of each comparison's own total
         avg_log10 = log10((A + B) / 2)) %>%
  ungroup() %>%
  mutate(direction = case_when(log2FC >  fc_cutoff ~ paste("Higher in", B_name),
                               log2FC < -fc_cutoff ~ paste("Higher in", A_name),
                               TRUE                ~ "Similar"),
         short     = str_trunc(Name, 40))

# ---- 3. Per-comparison figures --------------------------------------
for (i in seq_len(nrow(comparisons))) {

  cmp  <- comparisons$comparison[i]
  A_nm <- comparisons$A[i]
  B_nm <- comparisons$B[i]
  tag  <- tag_of(cmp)
  df   <- filter(all_calc, comparison == cmp)

  message(sprintf("[%s] %d compounds detected in both cultures.", cmp, nrow(df)))

  cols <- setNames(c(culture_cols[[B_nm]], culture_cols[[A_nm]], "grey70"),
                    c(paste("Higher in", B_nm), paste("Higher in", A_nm), "Similar"))

  top_lab <- df %>% slice_max(abs(log2FC), n = 8, with_ties = FALSE)

  # Reference lines follow the same normalization as log2FC (share of total
  # signal), so the dotted lines match the colouring exactly.
  offset <- log10(sum(df$B) / sum(df$A))

  # -- Scatter plot: A vs B --
  p_scatter <- ggplot(df, aes(A, B, color = direction)) +
    geom_abline(slope = 1, intercept = offset, linetype = "dashed") +
    geom_abline(slope = 1, intercept = offset + fc_cutoff * log10(2), linetype = "dotted") +
    geom_abline(slope = 1, intercept = offset - fc_cutoff * log10(2), linetype = "dotted") +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(data = top_lab, aes(label = short), color = "black",
                    size = 3, max.overlaps = Inf, min.segment.length = 0, seed = 42) +
    scale_color_manual(values = cols) +
    scale_x_log10(labels = label_log()) +
    scale_y_log10(labels = label_log()) +
    labs(title = sprintf("%s vs %s intensities", B_nm, A_nm),
         subtitle = "dashed = equal share of total signal; dotted = 4-fold difference",
         x = sprintf("%s (OD-normalized intensity)", A_nm),
         y = sprintf("%s (OD-normalized intensity)", B_nm), color = NULL) +
    base_theme + theme(legend.position = "none")

  ggsave(sprintf("EPI300_scatter_%s.png", tag), p_scatter, width = 8, height = 7, dpi = 300)

  # -- MA plot --
  p_ma <- ggplot(df, aes(avg_log10, log2FC, color = direction)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(data = top_lab, aes(label = short), color = "black",
                    size = 3, max.overlaps = Inf, min.segment.length = 0, seed = 42) +
    scale_color_manual(values = cols) +
    labs(title = sprintf("MA plot: %s", cmp),
         subtitle = "fold change vs average abundance",
         x = expression(log[10]~average~intensity),
         y = bquote(log[2]~fold~change~"("*.(B_nm)*" / "*.(A_nm)*")")) +
    base_theme + theme(legend.position = "none")

  ggsave(sprintf("EPI300_MA_plot_%s.png", tag), p_ma, width = 9, height = 7, dpi = 300)

  # -- Ranked bar chart: top up / down --
  n_side <- 15
  bars <- bind_rows(slice_max(df, log2FC, n = n_side, with_ties = FALSE),
                    slice_min(df, log2FC, n = n_side, with_ties = FALSE)) %>%
    distinct(Name, .keep_all = TRUE) %>%
    mutate(short = fct_reorder(short, log2FC))

  p_bar <- ggplot(bars, aes(log2FC, short, fill = log2FC > 0)) +
    geom_col() +
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = setNames(c(culture_cols[[B_nm]], culture_cols[[A_nm]]),
                                        c(TRUE, FALSE))) +
    labs(title = sprintf("Top %d up / down compounds: %s", n_side, cmp),
         x = bquote(log[2]~fold~change~"("*.(B_nm)*" / "*.(A_nm)*")"), y = NULL) +
    base_theme + theme(legend.position = "none",
                       axis.text.y = element_text(size = 9))

  ggsave(sprintf("EPI300_top_up_down_bars_%s.png", tag), p_bar, width = 10, height = 9, dpi = 300)

  # -- Heatmap: top compounds by |log2FC| --
  n_heat <- 30
  heat <- df %>%
    slice_max(abs(log2FC), n = n_heat, with_ties = FALSE) %>%
    mutate(short = fct_reorder(short, log2FC)) %>%
    pivot_longer(c(A, B), names_to = "slot", values_to = "intensity") %>%
    mutate(Culture = if_else(slot == "A", A_nm, B_nm),
           Culture = factor(Culture, levels = c(A_nm, B_nm)))

  p_heat <- ggplot(heat, aes(Culture, short, fill = log10(intensity))) +
    geom_tile(color = "white") +
    geom_text(aes(label = comma(round(intensity)),
                  color = log10(intensity) > 5), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                         name = expression(log[10]~intensity)) +
    labs(title = sprintf("Heatmap: top %d compounds by fold change (%s)", n_heat, cmp),
         subtitle = sprintf("ordered by log2FC (%s / %s); numbers = intensity", B_nm, A_nm),
         x = NULL, y = NULL) +
    base_theme + theme(panel.grid = element_blank(),
                       axis.text.y = element_text(size = 9))

  ggsave(sprintf("EPI300_heatmap_%s.png", tag), p_heat, width = 9, height = 9, dpi = 300)

  # -- Most abundant compounds, both cultures --
  n_abun <- 15
  abun <- df %>%
    slice_max(A + B, n = n_abun, with_ties = FALSE) %>%
    mutate(short = fct_reorder(short, A + B)) %>%
    pivot_longer(c(A, B), names_to = "slot", values_to = "intensity") %>%
    mutate(Culture = if_else(slot == "A", A_nm, B_nm),
           Culture = factor(Culture, levels = c(A_nm, B_nm)))

  p_abun <- ggplot(abun, aes(intensity, short, fill = Culture)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = culture_cols) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(title = sprintf("Top %d most abundant compounds: %s", n_abun, cmp),
         x = "OD-normalized intensity", y = NULL, fill = NULL) +
    base_theme + theme(legend.position = "bottom",
                       axis.text.y = element_text(size = 9))

  ggsave(sprintf("EPI300_top_abundant_%s.png", tag), p_abun, width = 10, height = 8, dpi = 300)
}

# ---- 4. Combined overview: scatter & MA, all comparisons side by side --
# Build one global colour map so "Higher in Ara" etc. is the same colour
# in every facet (keyed off culture_cols), instead of ggplot's default
# per-plot discrete palette.
all_directions <- unique(all_calc$direction)
dir_cols <- setNames(
  ifelse(all_directions == "Similar", "grey70",
        culture_cols[str_remove(all_directions, "^Higher in ")]),
  all_directions
)

p_scatter_all <- ggplot(all_calc, aes(A, B, color = direction)) +
  geom_point(size = 1.6, alpha = 0.75) +
  scale_x_log10(labels = label_log()) +
  scale_y_log10(labels = label_log()) +
  scale_color_manual(values = dir_cols) +
  facet_wrap(~ comparison, nrow = 1, scales = "free") +
  labs(title = "A vs B intensities, all comparisons",
       subtitle = "x = first-named condition, y = second-named condition (see facet title)",
       x = "Condition A (OD-normalized intensity)",
       y = "Condition B (OD-normalized intensity)", color = NULL) +
  base_theme + theme(legend.position = "bottom",
                     strip.text = element_text(face = "bold"))

ggsave("EPI300_scatter_all.png", p_scatter_all, width = 17, height = 5.5, dpi = 300)

p_ma_all <- ggplot(all_calc, aes(avg_log10, log2FC, color = direction)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
  geom_point(size = 1.6, alpha = 0.75) +
  scale_color_manual(values = dir_cols) +
  facet_wrap(~ comparison, nrow = 1, scales = "free") +
  labs(title = "MA plots, all comparisons",
       subtitle = "fold change vs average abundance (B / A per facet title)",
       x = expression(log[10]~average~intensity),
       y = expression(log[2]~fold~change), color = NULL) +
  base_theme + theme(legend.position = "bottom",
                     strip.text = element_text(face = "bold"))

ggsave("EPI300_MA_all.png", p_ma_all, width = 17, height = 5.5, dpi = 300)

message("Saved per-comparison scatter, MA, top up/down bars, heatmap, top abundant ",
        "(x4 comparisons), plus EPI300_scatter_all.png and EPI300_MA_all.png.")
