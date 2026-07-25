#### 05_Prevalence_Analysis.R — Weighted disease prevalence + comparison plots
#### Run inside the AoU Researcher Workbench, AFTER 04_Global_RAILS.R has
#### produced the per-participant weights (global_rails_weights.csv).
####
#### Estimates GBD-category prevalence under each weighting method, compares
#### to US reference prevalence, and produces the national and by-region
#### figures. Assumes `query()` (BigQuery helper) and `codes` (ICD lookup
#### table with columns Disease / icd9 / icd10) are available in the session.

library(tidyverse)
library(survey)

########################################################################
## Load per-participant weights
##   Columns: person_id, covariates, region, and the weight ladder
##   w_unweighted, w_cal1, w_cal2, w_nps1, w_nps2, w_nps1_rake,
##   w_nps2_rake, w_rails
########################################################################

my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
system(paste0("gsutil cp ", my_bucket, "/data/global_rails_weights.csv ."), intern = TRUE)
raking.wts <- read_csv("global_rails_weights.csv")

########################################################################
## Flag each GBD disease category from ICD9/ICD10 condition codes
##   A participant is a case if they have the condition on >1 distinct date.
##   (Generalizes the original 15 copy-pasted blocks into one loop.)
########################################################################

## Build the BETWEEN-range SQL clause for one disease's icd9 + icd10 ranges
build_icd_search <- function(disease_row) {
  parse_ranges <- function(codes_str) {
    do.call(rbind, strsplit(codes_str, split = ", ")) %>%
      t() %>% as.data.frame() %>% setNames("code") %>%
      tidyr::separate(code, into = c("first", "last"), sep = "-", fill = "right") %>%
      mutate(last = ifelse(is.na(last), first, last)) %>%
      distinct()
  }
  ranges <- rbind(parse_ranges(disease_row$icd9), parse_ranges(disease_row$icd10))
  clauses <- sprintf("(condition_source_value BETWEEN '%s' AND '%s')",
                     ranges$first, ranges$last)
  paste(clauses, collapse = " OR ")
}

## Diseases to flag: internal name -> label in the `codes` table
disease_defs <- tibble::tribble(
  ~var,        ~disease_label,
  "digestive", "Digestive diseases",
  "neuro",     "Neurological disorders",
  "cvd",       "Cardiovascular diseases",
  "skin",      "Skin and subcutaneous diseases",
  "cancer",    "Neoplasms",
  "preg",      "Maternal and neonatal disorders",
  "mental",    "Mental disorders",
  "mskel",     "Musculoskeletal disorders",
  "diab",      "Diabetes mellitus",
  "resp",      "Respiratory infections and tuberculosis",
  "nutr",      "Nutritional deficiencies",
  "kidney",    "Chronic kidney disease",
  "vis",       "Blindness and vision loss",
  "abuse",     "Substance use disorders",
  "hiv",       "hiv"
)

for (i in seq_len(nrow(disease_defs))) {
  var   <- disease_defs$var[i]
  label <- disease_defs$disease_label[i]
  search <- build_icd_search(filter(codes, Disease == label))

  q <- str_glue("
    SELECT person_id AS Person,
           COUNT(DISTINCT condition_start_date) AS N
    FROM condition_occurrence AS co
    JOIN concept AS mc ON (condition_source_concept_id = mc.concept_id)
    WHERE mc.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    AND ({search})
    GROUP BY Person
  ")

  cases <- query(q) %>% filter(N > 1)
  raking.wts[[var]] <- ifelse(raking.wts$person_id %in% cases$Person, 1, 0)
  message(label, ": ", sum(raking.wts[[var]]), " cases")
}

########################################################################
## Weighted prevalence with a design-based CI
########################################################################

wt.prev <- function(w, y) {
  ok  <- !is.na(w) & !is.na(y)
  w   <- w[ok]; y <- y[ok]
  wm   <- sum(w * y) / sum(w)
  wvar <- sum(w^2 * (y - wm)^2) / (sum(w)^2)
  c(EST = wm,
    LB  = wm - qnorm(0.975) * sqrt(wvar),
    UB  = wm + qnorm(0.975) * sqrt(wvar))
}

## Method ladder: internal name (used downstream) -> weight column
method_cols <- c(
  raw              = "w_unweighted",
  svy_oneway_rake  = "w_cal1",
  ps_oneway        = "w_nps1",
  ps_oneway_rake   = "w_nps1_rake",
  ps_twoway        = "w_nps2",
  ps_twoway_rake   = "w_nps2_rake",
  ps_threeway_rake = "w_rails"
)

## Prevalence (per 1e6) for one disease under every method, wide with LB/UB
calc.wt.prev <- function(df, dis) {
  parts <- lapply(names(method_cols), function(m) {
    est <- wt.prev(df[[method_cols[m]]], df[[dis]]) * 1e6
    setNames(as.list(est), paste0(m, c("", "_LB", "_UB")))
  })
  cbind(data.frame(dis = dis), as.data.frame(do.call(c, parts)))
}

########################################################################
## National prevalence table
########################################################################

all <- bind_rows(lapply(disease_defs$var, function(d) calc.wt.prev(raking.wts, d)))

## US reference prevalence (per 1e6), aligned to disease_defs order
us.prev <- c(0.337497, 0.50263, 0.160279, 0.29664, 0.104958, 0.013393,
             0.176802, 0.49361, 0.156573, 0.17518, 0.055389, 0.161258,
             0.040259, 0.067958, 0.265634) * 1e6
all <- cbind(all, us.prev)

########################################################################
## Reshape to long: one row per disease x method with EST/Lower/Upper
########################################################################

dis_levels <- disease_defs$var
dis_labels <- c("Digestive", "Neurological", "Cardiovascular", "Skin", "Cancer",
                "Maternal/Neonatal", "Mental", "Musculoskeletal", "Diabetes",
                "Respiratory/TB", "Nutritional", "Kidney", "Vision",
                "Substance Abuse", "HIV/STIs")

method_levels <- c("us.prev", names(method_cols))
method_labels <- c("US Prevalence", "Raw", "Survey Oneway Rake", "PS Oneway",
                   "PS Oneway Rake", "PS Twoway", "PS Twoway Rake",
                   "PS Var. Sel. \n Threeway Stepwise \n Rake")

## Final display labels used in the figures
final_levels <- c("US Prevalence", "Raw", "Survey Oneway Rake", "PS Oneway",
                  "PS Oneway Rake", "PS Twoway", "PS Twoway Rake",
                  "PS Var. Sel. \n Threeway Stepwise \n Rake")
final_labels <- c("US", "naive", "cal-1", "nps-1", "nps-cal-1",
                  "nps-2", "nps-cal-2", "RAILS")

## Helper: pivot the point estimates / a bound column set to long form
to_long <- function(df, cols, value_name, lvls, lbls, region_name = NULL) {
  out <- df %>%
    select(dis, all_of(cols)) %>%
    tidyr::pivot_longer(-dis, names_to = "method", values_to = value_name) %>%
    mutate(
      method = factor(method, levels = lvls, labels = lbls),
      dis    = factor(dis,    levels = dis_levels, labels = dis_labels)
    )
  if (!is.null(region_name)) out$region <- region_name
  out
}

build_long <- function(df) {
  est <- to_long(df, c("us.prev", names(method_cols)),
                 "prevalence", method_levels, method_labels)
  lb  <- to_long(df, paste0(names(method_cols), "_LB"),
                 "Lower", paste0(names(method_cols), "_LB"), method_labels[-1])
  ub  <- to_long(df, paste0(names(method_cols), "_UB"),
                 "Upper", paste0(names(method_cols), "_UB"), method_labels[-1])
  est %>%
    left_join(lb, by = c("dis", "method")) %>%
    left_join(ub, by = c("dis", "method"))
}

## Disease order used in the figures (by overall burden / grouping)
fig_dis_order <- c("Maternal/Neonatal", "Vision", "Nutritional", "Substance Abuse",
                   "Cancer", "Diabetes", "Cardiovascular", "Kidney",
                   "Respiratory/TB", "Mental", "HIV/STIs", "Skin",
                   "Digestive", "Musculoskeletal", "Neurological")

all.long <- build_long(all) %>%
  mutate(method = factor(method, levels = final_levels, labels = final_labels),
         dis    = factor(dis, levels = fig_dis_order)) %>%
  arrange(dis, desc(method)) %>%
  mutate(ordering = factor(seq_len(n()))) %>%
  group_by(dis) %>%
  mutate(prev.us.ratio = prevalence / prevalence[method == "US"],
         LB.us.ratio    = Lower      / prevalence[method == "US"],
         UB.us.ratio    = Upper      / prevalence[method == "US"]) %>%
  ungroup()

########################################################################
## National figures
########################################################################

options(scipen = 1e6)

## (1) All methods, absolute prevalence.
## Distinct shape per method (mix of filled 22-25 and open 0/1/3/4 markers)
## so methods are distinguishable in grayscale / colorblind print, not all
## filled circles. color is mapped alongside fill so the open shapes — which
## ignore fill — still carry the method color.
method_cols_fig <- c("red", "green", "turquoise1", "dodgerblue1",
                     "#B600FF", "#FF00DB", "chocolate", "orange")
method_shapes   <- c(25, 24, 0, 1, 3, 23, 4, 22)   # US, naive, cal-1, nps-1,
                                                   # nps-cal-1, nps-2, nps-cal-2, RAILS
method_sizes    <- c(1.4, 1.4, 1, 1, 1, 1, 1, 1)

## Point size is set per method OUTSIDE aes() via a lookup, so `size` is not
## a mapped aesthetic and produces no legend of its own. The three mapped
## scales (fill/colour/shape) share the name "Method" and identical default
## guides, so ggplot merges them into ONE legend whose keys carry each
## series' real shape and colour — no duplicated colour keys, no
## override.aes warnings.
p_all <- ggplot(all.long, aes(x = prevalence, y = ordering)) +
  geom_point(aes(fill = method, color = method, shape = method),
             size = method_sizes[as.integer(all.long$method)]) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 1, linewidth = 0.25) +
  scale_fill_manual("Method",  values = method_cols_fig) +
  scale_color_manual("Method", values = method_cols_fig) +
  scale_shape_manual("Method", values = method_shapes) +
  labs(x = "Prevalence (per 1 million)", y = "Disease Category") +
  theme_bw() +
  theme(axis.text.y = element_blank(), legend.position = "top",
        legend.key.size = unit(0.6, "cm")) +
  facet_grid(dis ~ ., scales = "free", space = "free_y", switch = "y") +
  theme(strip.placement = "outside", panel.spacing = unit(0, "in"),
        strip.background.y = element_rect(fill = "white", color = "gray75"),
        strip.text.y.left = element_text(angle = 0))
print(p_all)

## Publication-quality export of Figure 1.
## 600 dpi TIFF with LZW compression is the usual journal requirement; the
## PDF is vector (infinitely scalable) and is preferred when accepted.
## Sized for a full-page single-column figure — adjust width/height if the
## journal specifies exact dimensions.
ggsave("Fig1_prevalence_all_methods.tiff", plot = p_all,
       width = 7.5, height = 10, units = "in", dpi = 600,
       device = "tiff", compression = "lzw")

ggsave("Fig1_prevalence_all_methods.pdf", plot = p_all,
       width = 7.5, height = 10, units = "in", device = cairo_pdf)

## High-res PNG. Blurry text in the manuscript is almost always type rendered
## too small for the pixel grid, not a dpi shortfall — so bump BOTH the
## resolution (to 1200 dpi) and the base font size, and render text crisply
## with the ragg PNG device (falls back to the default device if ragg is
## absent). 7.5 x 10 in at 1200 dpi = 9000 x 12000 px.
p_all_print <- p_all + theme(
  text            = element_text(size = 11),
  axis.text.x     = element_text(size = 9),
  legend.text     = element_text(size = 10),
  legend.title    = element_text(size = 11),
  strip.text.y.left = element_text(angle = 0, size = 9)
)

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave("Fig1_prevalence_all_methods.png", plot = p_all_print,
         width = 7.5, height = 10, units = "in", dpi = 1200, device = ragg::agg_png)
} else {
  ggsave("Fig1_prevalence_all_methods.png", plot = p_all_print,
         width = 7.5, height = 10, units = "in", dpi = 1200, type = "cairo")
}

# system(paste0("gsutil cp ./Fig1_prevalence_all_methods.* ", my_bucket, "/data/"), intern = TRUE)

## (2) US vs naive vs RAILS only
p_key <- all.long %>%
  filter(method %in% c("US", "naive", "RAILS")) %>%
  mutate(method = factor(method, levels = c("US", "naive", "RAILS"),
                         labels = c("US", "Unweighted", "RAILS"))) %>%
  ggplot(aes(x = prevalence, y = ordering)) +
  geom_point(aes(fill = method, shape = method, size = method)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 1, linewidth = 0.25) +
  scale_fill_manual(values = c("red", "green", "orange")) +
  scale_shape_manual(values = c(25, 23, 21)) +
  scale_size_manual(values = c(2, 2, 1.5)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = "none", size = "none") +
  labs(x = "Prevalence (per 1 million)", y = "Disease Category", fill = "Method") +
  theme_bw() +
  theme(axis.text.y = element_blank(), legend.position = "top") +
  facet_grid(dis ~ ., scales = "free", space = "free_y", switch = "y") +
  theme(strip.placement = "outside", panel.spacing = unit(0, "in"),
        strip.background.y = element_rect(fill = "white", color = "gray75"),
        strip.text.y.left = element_text(angle = 0))
print(p_key)

## (3) Ratio to US prevalence (all methods except US)
p_ratio <- all.long %>%
  filter(method != "US") %>% droplevels() %>%
  ggplot(aes(x = prev.us.ratio, y = ordering)) +
  geom_point(aes(fill = method, shape = method, size = method)) +
  geom_errorbar(aes(xmin = LB.us.ratio, xmax = UB.us.ratio), width = 1, linewidth = 0.25) +
  geom_vline(xintercept = 1, color = "gray20", linetype = "dashed") +
  scale_fill_manual(values = c("green", "turquoise1", "dodgerblue1",
                               "#B600FF", "#FF00DB", "chocolate", "orange")) +
  scale_shape_manual(values = c(23, 21, 21, 21, 21, 21, 21)) +
  scale_size_manual(values = c(2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = "none", size = "none") +
  labs(x = "Ratio to US Prevalence", y = "Disease Category", fill = "Method") +
  theme_bw() +
  theme(axis.text.y = element_blank(), legend.position = "top") +
  facet_grid(dis ~ ., scales = "free", space = "free_y", switch = "y") +
  theme(strip.placement = "outside", panel.spacing = unit(0, "in"),
        strip.background.y = element_rect(fill = "white", color = "gray75"),
        strip.text.y.left = element_text(angle = 0))
print(p_ratio)

########################################################################
## By-region prevalence
##   region already lives in raking.wts (from 04_Global_RAILS.R); no re-query.
########################################################################

calc.wt.prev.region <- function(df, dis, reg) {
  calc.wt.prev(filter(df, region == reg), dis)
}

region_levels <- c("South", "Midwest", "Northeast", "West")

full <- bind_rows(lapply(region_levels, function(reg) {
  tab <- bind_rows(lapply(disease_defs$var,
                          function(d) calc.wt.prev.region(raking.wts, d, reg)))
  tab <- cbind(tab, us.prev)
  build_long(tab) %>% mutate(region = reg)
})) %>%
  mutate(method = factor(method, levels = final_levels, labels = final_labels),
         dis    = factor(dis, levels = fig_dis_order),
         region = factor(region, levels = region_levels))

## full2: methods on the y-axis within each region facet
full2 <- full %>%
  group_by(dis, region) %>%
  mutate(prev.us.ratio = prevalence / prevalence[method == "US"],
         LB.us.ratio    = Lower      / prevalence[method == "US"],
         UB.us.ratio    = Upper      / prevalence[method == "US"]) %>%
  ungroup() %>%
  filter(!(method %in% c("nps-1", "nps-cal-1")), method != "US") %>%
  group_by(region) %>%
  arrange(desc(method)) %>%
  mutate(ordering = factor(seq_len(n()))) %>%
  ungroup() %>%
  droplevels()

p_region_method <- ggplot(full2, aes(x = prevalence, y = ordering)) +
  geom_point(aes(fill = method, shape = method, size = method)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 1, linewidth = 0.15) +
  scale_fill_manual(values = c("green", "turquoise1", "#FF00DB", "chocolate", "orange")) +
  scale_shape_manual(values = c(23, 21, 21, 21, 21)) +
  scale_size_manual(values = c(1.5, 1, 1, 1, 1)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2)),
         shape = "none", size = "none") +
  labs(x = "Prevalence (per 1 million)", y = "Disease Category", fill = "Method") +
  lims(x = c(-28000, 700000)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 4), legend.text = element_text(size = 6),
        axis.title = element_text(size = 8), legend.title = element_text(size = 8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top") +
  facet_grid(dis ~ region, scales = "free", space = "fixed", switch = "y") +
  theme(strip.placement = "outside", panel.spacing = unit(0, "in"),
        strip.background.y = element_rect(fill = "white", color = "gray75"),
        strip.text.y.left = element_text(angle = 0, size = 6),
        strip.text.x = element_text(size = 6))
print(p_region_method)

## full3: regions colored within each method facet
full3 <- full %>%
  group_by(dis, region) %>%
  mutate(prev.us.ratio = prevalence / prevalence[method == "US"],
         LB.us.ratio    = Lower      / prevalence[method == "US"],
         UB.us.ratio    = Upper      / prevalence[method == "US"]) %>%
  ungroup() %>%
  filter(!(method %in% c("nps-1", "nps-cal-1")), method != "US") %>%
  group_by(method) %>%
  mutate(ordering = factor(seq_len(n()))) %>%
  ungroup() %>%
  droplevels()

p_method_region <- ggplot(full3, aes(x = prevalence, y = ordering)) +
  geom_point(aes(fill = region, shape = region, size = region)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 1, linewidth = 0.15) +
  scale_fill_manual(values = c("red", "green", "blue", "orange")) +
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_size_manual(values = c(1, 1, 1, 1)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2)),
         shape = "none", size = "none") +
  labs(x = "Prevalence (per 1 million)", y = "Disease Category", fill = "Region") +
  lims(x = c(-28000, 700000)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 4), legend.text = element_text(size = 6),
        axis.title = element_text(size = 8), legend.title = element_text(size = 8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top") +
  facet_grid(dis ~ method, scales = "free", space = "fixed", switch = "y") +
  theme(strip.placement = "outside", panel.spacing = unit(0, "in"),
        strip.background.y = element_rect(fill = "white", color = "gray75"),
        strip.text.y.left = element_text(angle = 0, size = 6),
        strip.text.x = element_text(size = 6))
print(p_method_region)

########################################################################
## Save the prevalence table
########################################################################

write_excel_csv(all,  "prevalence_national.csv")
write_excel_csv(full, "prevalence_by_region.csv")
system(paste0("gsutil cp ./prevalence_national.csv ",  my_bucket, "/data/"), intern = TRUE)
system(paste0("gsutil cp ./prevalence_by_region.csv ", my_bucket, "/data/"), intern = TRUE)
