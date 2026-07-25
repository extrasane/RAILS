#### 06_Phecode_Categories.R — Phecode-grouped disease prevalence figures
#### Run inside the AoU Researcher Workbench.
####
#### Maps AoU disease categories to phecodes and phecode groups, then plots
#### AoU vs US prevalence for every disease, colored by phecode grouping.
####
#### Inputs (workspace bucket):
####   raking_wts_w_diseases_*.csv  — person IDs + weights + disease indicators
####   Disease_Categories.csv       — disease -> ICD
####   phecodeX_unrolled_ICD_CM.csv — ICD9/10 -> phecode
####   PhecodeX_categories.csv      — phecode/ICD -> category name
####   AoU_vs_US.csv                — US reference prevalence per disease
####
#### Weight columns are set by WT_UNW / WT_RAILS below.

library(tidyverse)
library(ggplot2)
library(ggrepel)

`%nin%` <- Negate(`%in%`)

## Weight columns: unweighted benchmark and RAILS
WT_UNW   <- "w_unweighted"
WT_RAILS <- "w_rails"

########################################################################
## Load input files from the workspace bucket
########################################################################

my_bucket   <- Sys.getenv("WORKSPACE_BUCKET")
bucket_path <- paste0(my_bucket, "/data/")

fetch <- function(fname, path = bucket_path) {
  system(paste0("gsutil cp ", path, fname, " ."), intern = TRUE)
  read_csv(fname, show_col_types = FALSE)
}

raking.wts.w.diseases <- fetch("raking_wts_w_diseases_2_code_requirement.csv")
map_to_phecode <- fetch("Disease_Categories.csv")
phecodex       <- fetch("phecodeX_unrolled_ICD_CM.csv")
phecode_cats   <- fetch("PhecodeX_categories.csv")
prev           <- fetch("AoU_vs_US.csv")

########################################################################
## Map diseases -> phecode categories
##   Some ICD codes only match after stripping a trailing ".0"; a handful
##   need explicit overrides. Diseases matching several categories are
##   resolved by majority vote, then by the manual list below.
########################################################################

## Both joins are intentionally many-to-many: one ICD code can map to several
## phecodes, and one phecode covers several ICD codes. The duplication is
## resolved by the majority vote + manual list below.
map_to_phecode2 <- map_to_phecode %>%
  select(Disease, ICD) %>%
  left_join(phecode_cats, by = "ICD", relationship = "many-to-many") %>%
  mutate(ICD.v2 = ifelse(grepl(".0", ICD) & is.na(Phecode),
                         str_split(ICD, "\\.0", simplify = TRUE)[, 1],
                         NA)) %>%
  mutate(ICD.v2 = case_when(
    Disease == "Chronic kidney disease due to other and unspecified causes" ~ "753.1",
    Disease == "Ectopic pregnancy"                                          ~ "633.9",
    Disease %in% c("Osteoarthritis", "Osteoarthritis other")                ~ "715.9",
    TRUE                                                                    ~ ICD.v2)) %>%
  left_join(phecode_cats %>%
              select(ICD, Phecode, PhecodeCategory) %>%
              rename(phecode.v2 = Phecode, PhecodeCategoryV2 = PhecodeCategory),
            by = c("ICD.v2" = "ICD"), relationship = "many-to-many")

cats <- map_to_phecode2 %>%
  mutate(category = ifelse(!is.na(PhecodeCategoryV2), PhecodeCategoryV2, PhecodeCategory)) %>%
  mutate(category = ifelse(Disease %in% c("Maternal and neonatal disorders",
                                          "Maternal disorders",
                                          "Maternal obstructed labor and uterine rupture",
                                          "Other maternal disorders"),
                           "Pregnancy", category)) %>%
  filter(!is.na(category)) %>%
  group_by(Disease, category) %>%
  summarize(N = length(category), .groups = "drop") %>%
  group_by(Disease) %>%
  slice_max(N) %>%
  ungroup()

## Diseases still tied after the majority vote — resolved manually
manual_adj <- cats %>%
  count(Disease) %>%
  filter(n > 1) %>%
  pull(Disease)
message("Diseases needing manual category assignment: ", length(manual_adj))

cats <- cats %>%
  mutate(keep = case_when(
    Disease %nin% manual_adj ~ 1,
    Disease == "Acute hepatitis E"                                    & category == "Infections"      ~ 1,
    Disease == "Alcoholic cardiomyopathy"                             & category == "Cardiovascular"  ~ 1,
    Disease == "Benign and in situ cervical and uterine neoplasms"    & category == "Neoplasms"       ~ 1,
    Disease == "Chronic kidney disease due to diabetes mellitus type 1" & category == "Genitourinary" ~ 1,
    Disease == "Chronic kidney disease due to diabetes mellitus type 2" & category == "Genitourinary" ~ 1,
    Disease == "Chronic kidney disease due to hypertension"           & category == "Genitourinary"   ~ 1,
    Disease == "Congenital musculoskeletal and limb anomalies"        & category == "Congenital"      ~ 1,
    Disease == "Diabetes mellitus type 1"                             & category == "Endocrine/Metab" ~ 1,
    Disease == "Fungal skin diseases"                                 & category == "Dermatological"  ~ 1,
    Disease == "G6PD deficiency"                                      & category == "Blood/Immune"    ~ 1,
    Disease == "Hemolytic disease and other neonatal jaundice"        & category == "Neonatal"        ~ 1,
    Disease == "Maternal hypertensive disorders"                      & category == "Pregnancy"       ~ 1,
    Disease == "Maternal sepsis and other maternal infections"        & category == "Pregnancy"       ~ 1,
    Disease == "Multiple sclerosis"                                   & category == "Neurological"    ~ 1,
    Disease == "Neonatal preterm birth"                               & category == "Neonatal"        ~ 1,
    Disease == "Other chromosomal abnormalities"                      & category == "Genetic"         ~ 1,
    Disease == "Other chronic respiratory diseases"                   & category == "Respiratory"     ~ 1,
    Disease == "Other intestinal infectious diseases"                 & category == "Infections"      ~ 1,
    Disease == "Psoriasis"                                            & category == "Dermatological"  ~ 1,
    Disease == "Rheumatoid arthritis"                                 & category == "Muscloskeletal"  ~ 1,
    TRUE ~ 0)) %>%
  filter(keep == 1) %>%
  select(-keep)

message("Diseases mapped to a category: ", n_distinct(cats$Disease))

########################################################################
## Weighted prevalence per disease (unweighted vs RAILS), with CIs
########################################################################

## Disease indicator columns span from the first to the last disease name
dis_first <- "HIV/AIDS and sexually transmitted infections"
dis_last  <- "Executions and police conflict"

long_wts <- raking.wts.w.diseases %>%
  select(person_id, all_of(c(WT_UNW, WT_RAILS)), all_of(dis_first):all_of(dis_last)) %>%
  rename(w_unw = all_of(WT_UNW), w_rails = all_of(WT_RAILS)) %>%
  tidyr::pivot_longer(cols = -c(person_id, w_unw, w_rails),
                      names_to = "Disease", values_to = "presence")

## Diseases with <20 cases are too sparse to estimate — excluded downstream
l20 <- long_wts %>%
  group_by(Disease) %>%
  summarize(N = sum(presence), .groups = "drop") %>%
  filter(N < 20)

prevs <- long_wts %>%
  group_by(Disease) %>%
  summarize(
    AoU_unw  = sum(w_unw * presence) / sum(w_unw),
    unw_wvar = sum(w_unw^2 * (presence - AoU_unw)^2) / (sum(w_unw)^2),
    unw_LB   = AoU_unw - qnorm(0.975) * sqrt(unw_wvar),
    unw_UB   = AoU_unw + qnorm(0.975) * sqrt(unw_wvar),
    AoU_w    = sum(w_rails * presence) / sum(w_rails),
    w_wvar   = sum(w_rails^2 * (presence - AoU_w)^2) / (sum(w_rails)^2),
    w_LB     = AoU_w - qnorm(0.975) * sqrt(w_wvar),
    w_UB     = AoU_w + qnorm(0.975) * sqrt(w_wvar),
    .groups  = "drop"
  ) %>%
  select(-c(unw_wvar, w_wvar))

prev_ci <- prevs %>% filter(Disease %nin% l20$Disease)

########################################################################
## Per-disease estimate table (per 1,000 population), unweighted vs RAILS.
## Prints only the longtable BODY ROWS to the console — paste them between
## the \endlastfoot and \end{longtable} of the existing table skeleton.
########################################################################

## Selected BY RULE: every disease with a non-zero estimate under both
## weightings, alphabetical, to 3 decimals. (`prevs` has already dropped the
## <20-case diseases via `l20`.) This filter affects the printed table ONLY —
## the prevalence plots below use `data` / `data2` and are unchanged.
est_table <- prevs %>%
  filter(AoU_unw > 0, AoU_w > 0) %>%
  arrange(Disease) %>%
  mutate(across(c(AoU_unw, unw_LB, unw_UB, AoU_w, w_LB, w_UB), ~ .x * 1000))

message("Estimate table: ", nrow(est_table), " diseases with non-zero estimates")

cat(
  with(est_table,
       sprintf("%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
               Disease, AoU_unw, unw_LB, unw_UB, AoU_w, w_LB, w_UB)),
  sep = "\n"
)
cat("\n")

########################################################################
## Assemble the analysis table: category + US prevalence + AoU estimates
########################################################################

data <- cats %>%
  left_join(prev %>% select(Disease, US_prev), by = "Disease") %>%
  left_join(prevs %>% select(Disease, AoU_unw, AoU_w), by = "Disease") %>%
  filter(!is.na(US_prev), Disease %nin% l20$Disease) %>%
  mutate(
    Ratio_unw      = AoU_unw / US_prev,
    Ratio_w        = AoU_w   / US_prev,
    Ratio_w_to_unw = AoU_w   / AoU_unw,
    ## How much of the unweighted bias RAILS removes (>0 = moved toward US)
    rel_bias_reduc = (AoU_w - AoU_unw) / (AoU_unw - US_prev),
    ## Change in absolute bias (<0 = RAILS is closer to US)
    rel_bias_diff  = (abs(AoU_w - US_prev) - abs(AoU_unw - US_prev)) / abs(AoU_unw - US_prev)
  ) %>%
  select(Disease, category, US_prev, AoU_unw, AoU_w,
         Ratio_unw, Ratio_w, Ratio_w_to_unw, rel_bias_reduc, rel_bias_diff)

message("Diseases in the final table: ", nrow(data))

########################################################################
## Category factor, plot ordering, and annotation flags
########################################################################

cat_levels <- c("Cardiovascular", "Congenital", "Endocrine/Metab", "Genetic",
                "Gastrointestinal", "Genitourinary", "Blood/Immune", "Infections",
                "Muscloskeletal", "Neoplasms", "Neurological", "Sense organs",
                "Neonatal", "Pregnancy", "Mental", "Respiratory",
                "Dermatological", "Symptoms")
cat_labels <- c("Circulatory", "Congenital", "Endocrine", "Genetic",
                "Gastrointestinal", "Genitourinary", "Haemotopoietic", "Infections",
                "Musculoskeletal", "Neoplasms", "Neurologic", "Sense organs",
                "Perinatal", "Pregnancy", "Psychiatric", "Respiratory",
                "Skin", "Symptoms")

annotate_diseases <- c(
  "Cervical cancer", "Acute hepatitis C", "Ischemic heart disease",
  "Hypertensive heart disease", "Diabetes mellitus type 2", "Appendicitis",
  "Periodontal diseases", "G6PD deficiency", "Thalassemias",
  "Drug-susceptible tuberculosis", "Endometriosis",
  "Alzheimer's disease and other dementias", "Acute lymphoid leukemia",
  "Colon and rectum cancer", "Breast cancer", "Prostate cancer", "Low back pain",
  "Tracheal, bronchus, and lung cancer",
  "Neonatal encephalopathy due to birth asphyxia and trauma",
  "Headache disorders", "Maternal hypertensive disorders", "Maternal hemorrhage",
  "Multiple sclerosis", "Bipolar disorder", "Anxiety disorders")

data <- data %>%
  mutate(category = factor(category, levels = cat_levels, labels = cat_labels)) %>%
  arrange(category, Disease) %>%
  mutate(ordering    = seq_len(n()),
         is_annotate = ifelse(Disease %in% annotate_diseases, "yes", "no"),
         ## Also label anything RAILS pushed above the unweighted estimate
         annotate_ratio_g1 = ifelse(Ratio_w_to_unw > 1, "yes", is_annotate))

## x-axis: one tick per category, placed at its first disease
labels <- levels(data$category)
breaks <- data %>% group_by(category) %>% summarize(ordering = min(ordering)) %>% pull(ordering)

## Category palette (indices into rainbow(21) for maximum separation)
vec <- c(18, 4, 11, 19, 6, 12, 21, 7, 13, 1, 8, 14, 2, 9, 15, 3, 10, 16)
cat_palette <- rainbow(21)[vec]

########################################################################
## Shared theme for the large publication scatter figures
########################################################################

theme_big <- theme_bw() +
  theme(
    legend.text     = element_text(size = 24, face = "bold"),
    legend.title    = element_text(size = 24, face = "bold"),
    legend.position = c(0.18, 1),
    legend.direction = "horizontal",
    axis.line       = element_line(colour = "black"),
    axis.title.x    = element_text(size = 36),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 36),
    axis.title.y    = element_text(size = 36),
    axis.text.y     = element_text(hjust = 1, size = 36),
    panel.border    = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.margin     = unit(c(1, 1, 1.5, 1.2), "cm"),
    axis.ticks.length = unit(0.5, "cm")
  )

## Category scatter: y_var against disease ordering, sized by prevalence
plot_category_scatter <- function(df, y_var, y_lab, hline = 0,
                                  annot_col = "is_annotate", annot_size = 8,
                                  size_var = "AoU_unw", size_lab = "All of Us Prevalence") {
  ggplot(df, aes(x = ordering, y = .data[[y_var]],
                 size = .data[[size_var]], color = category)) +
    geom_point() +
    scale_colour_manual(values = cat_palette) +
    scale_size_continuous(range = c(4, 24)) +
    scale_x_continuous(labels = labels, breaks = breaks) +
    ## base subsetting here — .data[[...]] is a tidy-eval pronoun and is only
    ## valid inside aes()/dplyr verbs, not in the data argument
    geom_label_repel(data = df[df[[annot_col]] == "yes", , drop = FALSE],
                     aes(label = Disease), angle = 45, label.size = NA,
                     fill = "white", size = annot_size, box.padding = 2.5,
                     max.overlaps = Inf, color = "black") +
    geom_hline(yintercept = hline, linetype = "dashed", color = "gray", linewidth = 1) +
    labs(size = size_lab, x = "Disease Category", y = y_lab) +
    guides(color = "none") +
    theme_big
}

########################################################################
## Figure 3 — Ratio of RAILS-weighted to unweighted prevalence
########################################################################

print(summary(data$Ratio_w_to_unw))

## Diseases whose weighted estimate exceeds the unweighted one — these are
## labelled in the figure via the annotate_ratio_g1 flag (preset annotation
## list PLUS every disease with ratio > 1)
message("Diseases with weighted/unweighted ratio > 1: ", sum(data$Ratio_w_to_unw > 1, na.rm = TRUE))
print(data %>% filter(Ratio_w_to_unw > 1) %>%
        select(Disease, category, AoU_unw, AoU_w, Ratio_w_to_unw) %>%
        arrange(desc(Ratio_w_to_unw)))

p_fig3 <- plot_category_scatter(
  data, "Ratio_w_to_unw",
  "Prevalence Ratio (Weighted AoU vs Unweighted AoU)", hline = 1,
  annot_col = "annotate_ratio_g1")
print(p_fig3)
ggsave("Figure 3.jpg", plot = p_fig3, width = 48, height = 30, dpi = 300)

## Manuscript version: also label every disease with ratio > 1, larger labels
p_fig3_edit <- plot_category_scatter(
  data, "Ratio_w_to_unw",
  "Prevalence Ratio (RAILS-Weighted vs. Unweighted Estimates)", hline = 1,
  annot_col = "annotate_ratio_g1", annot_size = 12)
print(p_fig3_edit)
ggsave("Edited Figure.jpg", plot = p_fig3_edit, width = 48, height = 30, dpi = 300)

########################################################################
## Figure 4 — Relative bias reduction
########################################################################

print(summary(data$rel_bias_reduc))

p_fig4 <- plot_category_scatter(data, "rel_bias_reduc", "Relative Bias Reduction")
print(p_fig4)
ggsave("Figure 4.jpg", plot = p_fig4, width = 48, height = 30, dpi = 300)

## 4A: outliers dropped so the bulk of the points are legible
outliers <- c("Testicular cancer", "Diabetes mellitus type 2")
p_fig4a <- plot_category_scatter(
  data %>% filter(Disease %nin% outliers) %>%
    mutate(is_annotate = ifelse(Disease == "Diabetes mellitus type 2", "no", is_annotate)),
  "rel_bias_reduc", "Relative Bias Reduction")
print(p_fig4a)
ggsave("Figure 4A.jpg", plot = p_fig4a, width = 48, height = 30, dpi = 300)

########################################################################
## Figure 5 — Relative difference in bias
########################################################################

print(summary(data$rel_bias_diff))

p_fig5 <- plot_category_scatter(
  data %>% filter(Disease %nin% outliers) %>%
    mutate(is_annotate = ifelse(Disease == "Diabetes mellitus type 2", "no", is_annotate)),
  "rel_bias_diff", "Relative Difference in Bias")
print(p_fig5)
ggsave("Figure 5.jpg", plot = p_fig5, width = 48, height = 30, dpi = 300)

p_fig5a <- plot_category_scatter(data, "rel_bias_diff", "Relative Difference in Bias")
print(p_fig5a)
ggsave("Figure 5A.jpg", plot = p_fig5a, width = 48, height = 30, dpi = 300)

########################################################################
## Volcano-style plots — weighted vs unweighted ratio to US prevalence
##   Points below the diagonal moved closer to the US value under RAILS.
########################################################################

data <- data %>%
  mutate(dist      = abs(log(Ratio_w)),
         direction = ifelse(abs(log(Ratio_w)) < abs(log(Ratio_unw)), "Closer", "Further"),
         side      = ifelse(Ratio_unw < 1, "<1", ">1"))

finite_ratios <- data %>% filter(is.finite(Ratio_unw), is.finite(Ratio_w))

volcano_panel <- function(df, side_val, lims, tick_labels, title) {
  ggplot(df %>% filter(side == side_val),
         aes(x = log10(Ratio_unw), y = log10(Ratio_w),
             color = category, fill = category, size = AoU_w, shape = direction)) +
    geom_point(show.legend = TRUE) +
    geom_abline(intercept = 0, slope = 1) +
    scale_size_continuous(range = c(2, 0.25), limits = c(0, 0.6)) +
    scale_shape_manual(values = c(25, 24)) +
    scale_colour_manual(values = cat_palette, drop = FALSE) +
    scale_fill_manual(values = cat_palette, drop = FALSE) +
    scale_x_continuous(limits = lims, n.breaks = 4, labels = tick_labels, expand = c(0, 0)) +
    scale_y_continuous(limits = lims, n.breaks = 4, labels = tick_labels, expand = c(0, 0)) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2),
           size = guide_legend(order = 3), fill = "none") +
    labs(x = "Unweighted Ratio", y = "Weighted Ratio",
         color = "Disease \nCategory",
         shape = "Weighted AoU Proximity \nto US Prevalence",
         size  = "Weighted AoU Prevalence",
         title = title) +
    theme_bw() +
    theme(aspect.ratio = 1, legend.position = "top",
          legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5))
}

plot_l1 <- volcano_panel(finite_ratios, "<1", c(-2.5, 0.25),
                         c("0.001", "0.01", "0.1", "1"), "Unweighted Ratio <1")
plot_g1 <- volcano_panel(finite_ratios, ">1", c(-0.25, 3),
                         c("1", "10", "100", "1000"), "Unweighted Ratio >1")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  print((plot_l1 + plot_g1) + plot_layout(ncol = 2, guides = "collect") &
          theme(legend.position = "bottom", legend.box = "vertical"))
} else if (requireNamespace("ggpubr", quietly = TRUE)) {
  print(ggpubr::ggarrange(plot_l1, plot_g1, ncol = 2, common.legend = TRUE))
} else {
  print(plot_l1); print(plot_g1)
}

## Single-panel faceted version
p_volcano_facet <- ggplot(finite_ratios,
       aes(x = log10(Ratio_unw), y = log10(Ratio_w),
           color = category, fill = category, size = dist, shape = direction)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_size_continuous(range = c(2, 0.25)) +
  scale_shape_manual(values = c(25, 24)) +
  scale_colour_manual(values = cat_palette) +
  scale_fill_manual(values = cat_palette) +
  scale_x_continuous(limits = c(-3, 3), n.breaks = 7,
                     labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3, 3), n.breaks = 7,
                     labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                     expand = c(0, 0)) +
  guides(shape = "none", fill = "none") +
  labs(x = "Unweighted Ratio", y = "Weighted Ratio", color = "Disease \nCategory") +
  facet_wrap(~side, scales = "free",
             labeller = labeller(side = c("<1" = "Unweighted Ratio <1",
                                          ">1" = "Unweighted Ratio >1"))) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "top",
        legend.key.size = unit(0.25, "cm"))
print(p_volcano_facet)

########################################################################
## Figure 2B / 2C — RAILS prevalence ratio vs US, with confidence intervals
########################################################################

prevs_w_ci <- long_wts %>%
  group_by(Disease) %>%
  summarize(
    AoU  = sum(w_rails * presence) / sum(w_rails),
    wvar = sum(w_rails^2 * (presence - AoU)^2) / (sum(w_rails)^2),
    LB   = AoU - qnorm(0.975) * sqrt(wvar),
    UB   = AoU + qnorm(0.975) * sqrt(wvar),
    .groups = "drop"
  )

data2 <- data %>%
  left_join(prevs_w_ci, by = "Disease") %>%
  mutate(Ratio_lower = LB / US_prev,
         Ratio_upper = UB / US_prev,
         ## floor at 0.001 so the log scale stays finite
         Ratio_lower = ifelse(Ratio_lower <= 0, 0.001, Ratio_lower),
         ## HIGHLIGHT: weighted AoU prevalence exceeds the US reference
         above_one = is.finite(Ratio_w) & Ratio_w > 1)

## Diseases where AoU is over-represented relative to the US population
message("Diseases with weighted prevalence ratio > 1: ", sum(data2$above_one))
print(data2 %>% filter(above_one) %>%
        select(Disease, category, US_prev, AoU_w, Ratio_w) %>%
        arrange(desc(Ratio_w)))

plot_ratio_ci <- function(size_var, size_lab, highlight = TRUE) {
  df  <- data2 %>% filter(is.finite(Ratio_w))
  hi  <- df[df$above_one, , drop = FALSE]
  ## label the pre-set annotation list plus every highlighted disease
  lab <- df[df$is_annotate == "yes" | df$above_one, , drop = FALSE]

  p <- ggplot(df, aes(x = ordering, y = log10(Ratio_w),
                      size = .data[[size_var]], color = category)) +
    geom_point(alpha = if (highlight) 0.35 else 1)

  if (highlight) {
    ## Redraw the ratio > 1 points fully opaque with a black outline so they
    ## stand out against the de-emphasised rest of the series
    p <- p +
      geom_point(data = hi, shape = 21, stroke = 1.5,
                 colour = "black", aes(fill = category), show.legend = FALSE) +
      scale_fill_manual(values = cat_palette)
  }

  p +
    geom_errorbar(aes(ymin = log10(Ratio_lower), ymax = log10(Ratio_upper)),
                  width = 1, linewidth = 0.25, alpha = if (highlight) 0.35 else 1) +
    geom_errorbar(data = hi, aes(ymin = log10(Ratio_lower), ymax = log10(Ratio_upper)),
                  width = 1, linewidth = 0.4, colour = "black") +
    scale_colour_manual(values = cat_palette) +
    scale_size_continuous(range = c(4, 24)) +
    scale_x_continuous(labels = labels, breaks = breaks) +
    scale_y_continuous(limits = c(-3, 3), n.breaks = 7,
                       labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                       expand = c(0, 0)) +
    geom_label_repel(data = lab, aes(label = Disease),
                     angle = 45, label.size = NA, fill = "white", size = 8,
                     box.padding = 2.5, max.overlaps = Inf, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 1) +
    labs(size = size_lab, x = "Disease Category",
         y = "Prevalence Ratio (Weighted AoU vs U.S National)") +
    guides(color = "none", fill = "none") +
    theme_big
}

p_fig2b <- plot_ratio_ci("AoU_unw", "All of Us Prevalence")
print(p_fig2b)
ggsave("Figure 2B.jpg", plot = p_fig2b, width = 48, height = 30, dpi = 300)

p_fig2c <- plot_ratio_ci("US_prev", "US Prevalence")
print(p_fig2c)
ggsave("Figure 2C.jpg", plot = p_fig2c, width = 48, height = 30, dpi = 300)

########################################################################
## Save tables and figures back to the bucket
########################################################################

write_excel_csv(data,    "phecode_category_prevalence.csv")
write_excel_csv(prev_ci, "phecode_prevalence_with_ci.csv")
system(paste0("gsutil cp ./phecode_category_prevalence.csv ", my_bucket, "/data/"), intern = TRUE)
system(paste0("gsutil cp ./phecode_prevalence_with_ci.csv ",  my_bucket, "/data/"), intern = TRUE)
system(paste0("gsutil cp ./Figure*.jpg ",       my_bucket, "/data/"), intern = TRUE)
system(paste0("gsutil cp './Edited Figure.jpg' ", my_bucket, "/data/"), intern = TRUE)
