# ============================================================
# PCA template
#
# TWO MODES, chosen automatically:
#
# A) TRUE REPLICATES (preferred, if available)
#    Provide a CSV named REPLICATE_FILE below with: first column =
#    compound name, then ONE COLUMN PER SAMPLE. Replicate columns must
#    share a group name followed by a number, e.g. noAra_1, noAra_2,
#    noAra_3, Ara_1, Ara_2, Ara_3, iAAL_1 ... Each point in the PCA is
#    one sample; replicates of the same condition should cluster
#    together if the conditions really differ. 95% ellipses are drawn
#    for any group with >= 4 replicates.
#
# B) NO REPLICATES AVAILABLE (current EPI300_pCC1_ATF1 dataset)
#    We only have the 4 pairwise comparison CSVs, each with ONE value
#    per compound per condition (no repeated cultures). These are
#    merged into a single compound x condition table (noAra, Ara,
#    iAAL, Ara + iAAL), and PCA is run with each CONDITION as a single
#    "sample" (4 points total). This shows how the four conditions
#    separate in compound space, but -- with n = 1 per group -- it
#    CANNOT show within-condition (replicate) variability, so no
#    ellipses are drawn. Treat it as a descriptive ordination, not a
#    statistical test. If you later get true replicate data, drop it
#    in as REPLICATE_FILE and mode (A) takes over automatically.
# ============================================================

library(tidyverse)
library(ggrepel)

REPLICATE_FILE <- "EPI300_replicates.csv"   # <- true replicate-level file, if you have one

# The 4 pairwise comparison files used as a fallback when no replicate
# file is present (same layout as the other EPI300 scripts: column 1 =
# compound name, column 2 = condition A, column 3 = condition B).
comparisons <- tribble(
  ~comparison,             ~file,                                    ~A,      ~B,
  "noAra vs Ara",          "EPI300_pCC1_ATF1_noAraAraOD.csv",         "noAra", "Ara",
  "noAra vs iAAL",         "EPI300_pCC1_ATF1_noAraiAAlOD.csv",        "noAra", "iAAL",
  "noAra vs Ara + iAAL",   "EPI300_pCC1_ATF1_noAraAraiAAlOD.csv",     "noAra", "Ara + iAAL",
  "Ara vs Ara + iAAL",     "EPI300_pCC1_ATF1_AraAraiAAlOD.csv",       "Ara",   "Ara + iAAL"
)

# ---- 1a. Mode A: load true replicate-level file ---------------------
load_replicates <- function(file) {
  raw <- read_csv(file, show_col_types = FALSE)
  names(raw)[1] <- "Name"

  raw %>%
    filter(!is.na(Name), Name != "0") %>%
    mutate(Name = make.unique(Name, sep = " #")) %>%
    column_to_rownames("Name") %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
}

# ---- 1b. Mode B: merge the 4 pairwise comparison files ---------------
# A blank cell means "not detected" (-> 0). A cell that is non-blank but
# still fails to parse (e.g. an Excel "#VALUE!" error) is dropped with a
# message rather than silently zeroed. Repeated compound names within a
# file are summed (treated as one compound's total signal) so each
# compound gets exactly one row in the merged matrix.
read_condition_values <- function(comparison, file, A_name, B_name) {
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
    transmute(Name, A = replace_na(A_num, 0), B = replace_na(B_num, 0)) %>%
    pivot_longer(c(A, B), names_to = "slot", values_to = "value") %>%
    mutate(Condition = if_else(slot == "A", A_name, B_name)) %>%
    select(Name, Condition, value)
}

build_condition_matrix <- function(comparisons) {
  long <- pmap_dfr(comparisons, function(comparison, file, A, B)
                     read_condition_values(comparison, file, A, B))

  wide <- long %>%
    group_by(Name, Condition) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%  # collapse repeated peaks
    pivot_wider(names_from = Condition, values_from = value, values_fill = 0) %>%
    column_to_rownames("Name")

  # give the "Ara + iAAL" condition a plotting-safe sample name
  names(wide) <- str_replace_all(names(wide), " \\+ ", "_")
  as.matrix(wide)
}

# ---- 2. Pick a mode and build the matrix -----------------------------
if (file.exists(REPLICATE_FILE)) {
  message(sprintf("Found %s -- using true replicate-level PCA (mode A).", REPLICATE_FILE))
  mat <- load_replicates(REPLICATE_FILE)
  has_replicates <- TRUE
} else {
  message("No replicate file found -- building one sample per CONDITION from the ",
          "4 pairwise comparison CSVs (mode B). These are single measurements, not ",
          "biological replicates: no ellipses will be drawn, and this PCA shows how ",
          "the 4 conditions differ in compound space, not within-condition spread.")
  mat <- build_condition_matrix(comparisons)
  has_replicates <- FALSE
}

mat[is.na(mat)] <- 0

# ---- 3. Transform & filter -------------------------------------------
mat <- log2(mat + 1)                        # tame the huge dynamic range
mat <- mat[apply(mat, 1, var) > 0, ]        # drop compounds that never vary

# ---- 4. PCA (samples = rows, compounds = variables) -------------------
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

n_pc <- min(2, ncol(pca$x))
scores <- as_tibble(pca$x[, seq_len(n_pc), drop = FALSE], rownames = "Sample") %>%
  mutate(Group = str_remove(Sample, "[_\\-\\.]?\\d+$"))

# ---- 5. Plot ------------------------------------------------------------
n_per_group <- min(table(scores$Group))
draw_ellipse <- has_replicates && n_per_group >= 4

subtitle <- if (has_replicates) {
  NULL
} else {
  "1 sample per condition (no biological replicates) -- ordination only, no ellipses"
}

p <- ggplot(scores, aes(PC1, PC2, color = Group)) +
  { if (draw_ellipse) stat_ellipse(level = 0.95, linetype = "dotted") } +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 3.5, show.legend = FALSE,
                  seed = 42, min.segment.length = 0.2) +
  labs(title = "PCA of volatilome profiles",
       subtitle = subtitle,
       x = sprintf("PC1 (%s%%)", var_expl[1]),
       y = if (n_pc >= 2) sprintf("PC2 (%s%%)", var_expl[2]) else NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        panel.grid.minor = element_blank())

print(p)
ggsave("EPI300_PCA.png", p, width = 8, height = 6.5, dpi = 300)

# ---- 6. Which compounds drive PC1 / PC2? --------------------------------
loadings <- as_tibble(pca$rotation[, seq_len(n_pc), drop = FALSE], rownames = "Name") %>%
  arrange(desc(abs(PC1)))
write_csv(loadings, "EPI300_PCA_loadings.csv")

message(sprintf("Saved EPI300_PCA.png (%d samples, %d compounds used) and EPI300_PCA_loadings.csv.",
                nrow(scores), nrow(mat)))
