# ============================================================
# Venn diagram, box plot and violin plot
# EPI300_pCC1_ATF1 volatilome: noAra, Ara, iAAL, Ara + iAAL
#
# Each comparison file has: column 1 = compound name, column 2 = culture A,
# column 3 = culture B (OD-normalized intensities; 0 / blank = not detected).
# Columns are read BY POSITION, so header spelling doesn't matter.
# Outputs: EPI300_venn.png, EPI300_boxplot.png, EPI300_violin.png
# ============================================================

library(tidyverse)
library(scales)

# ---- 1. What to compare -------------------------------------------
comparisons <- tribble(
  ~comparison,             ~file,                                    ~A,      ~B,
  "noAra vs Ara",          "EPI300_pCC1_ATF1_noAraAraOD.csv",         "noAra", "Ara",
  "noAra vs iAAL",         "EPI300_pCC1_ATF1_noAraiAAlOD.csv",        "noAra", "iAAL",
  "noAra vs Ara + iAAL",   "EPI300_pCC1_ATF1_noAraAraiAAlOD.csv",     "noAra", "Ara + iAAL",
  "Ara vs Ara + iAAL",     "EPI300_pCC1_ATF1_AraAraiAAlOD.csv",       "Ara",   "Ara + iAAL"
)

culture_cols <- c("noAra" = "blue", "Ara" = "red",
                  "iAAL" = "darkgreen", "Ara + iAAL" = "purple")

# ---- 2. Read each file into one long table ------------------------
read_comparison <- function(comparison, file, A, B) {
  raw <- read_csv(file, col_types = cols(.default = col_character()),
                  show_col_types = FALSE, name_repair = "minimal")[, 1:3]
  names(raw) <- c("Name", "A", "B")

  d <- raw %>%
    filter(!is.na(Name), Name != "0") %>%            # empty / placeholder names
    mutate(A_num = suppressWarnings(as.numeric(A)),
           B_num = suppressWarnings(as.numeric(B)),
           # A blank cell means "not detected" (-> 0 later). A cell that is
           # non-blank but still fails to parse (e.g. an Excel "#VALUE!"
           # error) is a real problem and gets dropped, not silently zeroed.
           is_bad = (is.na(A_num) & !is.na(A) & A != "") |
                    (is.na(B_num) & !is.na(B) & B != ""))

  n_bad <- sum(d$is_bad)
  if (n_bad > 0)
    message(sprintf("[%s] dropped %d row(s) with non-numeric (non-blank) values.",
                    comparison, n_bad))

  d %>%
    filter(!is_bad) %>%
    mutate(A = replace_na(A_num, 0),
           B = replace_na(B_num, 0)) %>%
    mutate(Name = make.unique(Name, sep = " #")) %>% # repeated names = separate peaks
    pivot_longer(c(A, B), names_to = "slot", values_to = "intensity") %>%
    mutate(culture = if_else(slot == "A", A, B),
           comparison = comparison) %>%
    select(comparison, Name, culture, intensity)
}

all_long <- pmap_dfr(comparisons, read_comparison) %>%
  mutate(comparison = factor(comparison, levels = comparisons$comparison))

# ---- 3. Venn diagram (presence / absence) -------------------------
# A compound counts as "found" in a culture if its intensity is > 0.
# Repeated names are merged so each compound is counted once.
venn_counts <- all_long %>%
  mutate(Name = str_remove(Name, " #\\d+$")) %>%
  group_by(comparison, Name, culture) %>%
  summarise(found = any(intensity > 0), .groups = "drop") %>%
  group_by(comparison, Name) %>%
  summarise(set = paste(sort(culture[found]), collapse = "|"),
            .groups = "drop") %>%
  filter(set != "") %>%
  left_join(comparisons %>% select(comparison, A, B), by = "comparison") %>%
  rowwise() %>%
  mutate(region = case_when(
    set == paste(sort(c(A, B)), collapse = "|") ~ "both",
    set == A ~ "onlyA",
    set == B ~ "onlyB")) %>%
  ungroup() %>%
  count(comparison, region) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  left_join(comparisons %>% select(comparison, A, B), by = "comparison") %>%
  mutate(total = onlyA + both + onlyB,
         comparison = factor(comparison, levels = comparisons$comparison))

# Some source files only list compounds detected in BOTH cultures (no unique
# hits for either side) -- flag those automatically instead of assuming which
# comparison it is, and mark them with an asterisk + footnote.
full_overlap <- setNames(venn_counts$onlyA == 0 & venn_counts$onlyB == 0,
                          as.character(venn_counts$comparison))
comparison_labels <- setNames(
  ifelse(full_overlap, paste0(names(full_overlap), "*"), names(full_overlap)),
  names(full_overlap)
)
has_full_overlap <- any(full_overlap)
overlap_caption <- "* this file only lists compounds detected in both cultures, so it shows 100% overlap by construction"

circle <- function(cx, r = 1, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(t), y = r * sin(t))
}

venn_circles <- venn_counts %>%
  rowwise() %>%
  reframe(comparison,
          side = rep(c("A", "B"), each = 200),
          culture = rep(c(A, B), each = 200),
          bind_rows(circle(-0.65), circle(0.65)))

pct <- function(n, tot) sprintf("%s\n(%.0f%%)", comma(n), 100 * n / tot)

venn_labels <- venn_counts %>%
  transmute(comparison,
            lab_A = pct(onlyA, total), lab_both = pct(both, total),
            lab_B = pct(onlyB, total), A, B) %>%
  pivot_longer(c(lab_A, lab_both, lab_B), names_to = "pos", values_to = "label") %>%
  mutate(x = case_when(pos == "lab_A" ~ -1.15, pos == "lab_both" ~ 0, TRUE ~ 1.15))

venn_names <- venn_counts %>%
  transmute(comparison, `-1.05` = A, `1.05` = B) %>%
  pivot_longer(-comparison, names_to = "x", values_to = "culture") %>%
  mutate(x = as.numeric(x))

p_venn <- ggplot() +
  geom_polygon(data = venn_circles,
               aes(x, y, group = interaction(comparison, side), fill = culture),
               alpha = 0.35, color = "black") +
  geom_text(data = venn_labels, aes(x, 0, label = label), size = 4.2, lineheight = 0.95) +
  geom_text(data = venn_names, aes(x, 1.25, label = culture, color = culture),
            fontface = "bold", size = 4.6, show.legend = FALSE) +
  scale_fill_manual(values = culture_cols, guide = "none") +
  scale_color_manual(values = culture_cols) +
  facet_wrap(~ comparison, nrow = 1, labeller = as_labeller(comparison_labels)) +
  coord_fixed(xlim = c(-2.1, 2.1), ylim = c(-1.1, 1.45), expand = FALSE) +
  labs(title = "Compounds detected in each culture",
       subtitle = "numbers = compounds (% of all compounds detected in either culture)",
       caption = if (has_full_overlap) overlap_caption else NULL) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 14)),
        plot.caption = element_text(size = 9, color = "grey30", hjust = 0.5),
        strip.text = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))

ggsave("EPI300_venn.png", p_venn, width = 19, height = 4.6, dpi = 300, bg = "white")

# ---- 4. Box plot & violin plot ------------------------------------
# Only compounds actually detected (> 0) in a culture are shown; zeros
# can't be placed on a log axis. n = number of compounds plotted.
dist_df <- all_long %>%
  filter(intensity > 0) %>%
  group_by(comparison, culture) %>%
  mutate(n = n(),
         xlab = sprintf("%s\n(n = %d)", culture, n)) %>%
  ungroup() %>%
  mutate(xlab = fct_reorder(xlab, match(culture, names(culture_cols))))

star_note <- if (has_full_overlap) "* this file only lists compounds detected in both cultures" else NULL

medians <- dist_df %>%
  group_by(comparison, xlab, culture) %>%
  summarise(med = median(intensity), .groups = "drop")

dist_theme <- theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.caption = element_text(size = 9, color = "grey30"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

y_scale <- scale_y_log10(labels = label_log(), expand = expansion(mult = c(0.05, 0.14)))
ylab <- "OD-normalized intensity (log scale)"

# Box plot (with individual compounds as points)
p_box <- ggplot(dist_df, aes(xlab, intensity, fill = culture)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.35, color = "grey20") +
  geom_label(data = medians, aes(y = Inf, label = paste0("median\n", comma(round(med)))),
             vjust = 1.1, size = 3.4, fill = "white", alpha = 0.85,
             label.size = 0, lineheight = 0.95) +
  scale_fill_manual(values = culture_cols) +
  y_scale +
  facet_wrap(~ comparison, nrow = 1, scales = "free_x", labeller = as_labeller(comparison_labels)) +
  labs(title = "Intensity distribution per culture",
       subtitle = "box = median and quartiles; dots = individual compounds",
       x = NULL, y = ylab, caption = star_note) +
  dist_theme

ggsave("EPI300_boxplot.png", p_box, width = 17, height = 6.5, dpi = 300)

# Violin plot (with a small box inside)
p_violin <- ggplot(dist_df, aes(xlab, intensity, fill = culture)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE, color = "grey20") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9) +
  scale_fill_manual(values = culture_cols) +
  y_scale +
  facet_wrap(~ comparison, nrow = 1, scales = "free_x", labeller = as_labeller(comparison_labels)) +
  labs(title = "Intensity distribution per culture",
       subtitle = "violin = shape of the distribution; inner box = median and quartiles",
       x = NULL, y = ylab, caption = star_note) +
  dist_theme

ggsave("EPI300_violin.png", p_violin, width = 17, height = 6.5, dpi = 300)

message("Saved: EPI300_venn.png, EPI300_boxplot.png, EPI300_violin.png")
print(venn_counts %>% select(comparison, A, B, onlyA, both, onlyB, total))
