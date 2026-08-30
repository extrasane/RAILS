#### BRCA_Survival_KM.R — Breast-cancer survival by BRCA / family history,
#### comparing three weighting schemes on the FEMALE AoU cohort:
####   1. Unweighted
####   2. Global RAILS      (w_rails,    the LIFO stepwise-raking G-RAILS weight)
####   3. Subgroup RAILS    (w_subrails, RAILS calibrated within sex = Female)
####
#### This reads the ALREADY-PREPARED analytic table (Data_for_Kaplan_Meier.csv),
#### which carries the survival window (start_age/end_age/ehr_length), the
#### phenotype flags, the demographics, and the precomputed IPW weight `w`
#### (balances the FH+WGS+EHR sub-cohort back to female AoU). Because that prep
#### is done, NO BigQuery / CDR access is needed here — we only attach the two
#### RAILS weights and draw the curves.
####
#### Survival weights (each combined with the precomputed IPW `w`, exactly as the
#### original notebook built wt.us = weight_lifo_rails * w):
####   wt_global = w_rails    * w
####   wt_sub    = w_subrails * w
####
#### Inputs (../data/):
####   Data_for_Kaplan_Meier.csv     (person_id, start_age, end_age, ehr_length,
####                                  has_FH_WGS_EHR, has_breastcancerDX,
####                                  has_BRCA1, has_BRCA2, has_BreastCancer_FH_any,
####                                  Race_simplified, Age, ..., w)
####   global_rails_weights.csv      (person_id, w_rails)
####   dt_sub_aou_femaleonly.csv     (person_id, w_subrails)   [by-sex sub-RAILS]
#### Outputs (../graph/):
####   km_three_schemes_row.png       KM panels, one per scheme (strata)
####   cuminc_three_schemes_row.png   cumulative-risk panels (1 - S), same strata
####   km_overlay_schemes.png         pooled KM overlay
####   cuminc_overlay_schemes.png     pooled cumulative-risk overlay
####   km_tail_trim_compare.png       tail, untrimmed vs winsorized weights
#### plus ../data/weighted_data_brca_bysex.csv

## Ensure a writable personal library, install any missing packages there.
## (In a conda-forge R env these are already present, so this is a no-op.)
user_lib <- path.expand(Sys.getenv("R_LIBS_USER", unset = file.path("~", "R", "library")))
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
for (p in c("tidyverse", "survival", "survminer"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, lib = user_lib, dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(survminer)
})

########################################################################
## Settings
########################################################################

KM_FILE <- "../data/Data_for_Kaplan_Meier.csv"   # alt: Cohort_C_02032026.csv (same schema + ancestry_recat)

########################################################################
## Load the prepared analytic table + derive has_BRCA (1 = BRCA1 or BRCA2)
########################################################################

km <- read_csv(KM_FILE, show_col_types = FALSE) %>%
  mutate(
    has_BRCA = as.integer(coalesce(has_BRCA1, 0) == 1 | coalesce(has_BRCA2, 0) == 1),
    ## has_FH_WGS_EHR may be stored as 0/1 or "No"/"Yes" — normalise to a flag
    fh_wgs_ehr = as.character(has_FH_WGS_EHR) %in% c("1", "Yes", "TRUE")
  )

message("Rows in prepared table: ", nrow(km),
        " | with all three sources: ", sum(km$fh_wgs_ehr, na.rm = TRUE))

########################################################################
## Attach RAILS weights: GLOBAL (w_rails) + SUBGROUP-by-sex (w_subrails),
## joined by person_id. Only two columns are read from the large weight files.
########################################################################

global_wts <- read_csv("../data/global_rails_weights.csv",
                       col_select = c(person_id, w_rails), show_col_types = FALSE) %>%
  distinct(person_id, .keep_all = TRUE)

sub_wts <- read_csv("../data/dt_sub_aou_femaleonly.csv",
                    col_select = c(person_id, w_subrails), show_col_types = FALSE) %>%
  distinct(person_id, .keep_all = TRUE)

km <- km %>%
  left_join(global_wts, by = "person_id") %>%
  left_join(sub_wts,    by = "person_id")

message("Matched global RAILS weight: ", sum(!is.na(km$w_rails)),
        " | subgroup RAILS weight: ",    sum(!is.na(km$w_subrails)))

########################################################################
## Analytic cohort + combined survival weights.
## Keep only participants with BOTH RAILS weights so all three curves are
## computed on the same sample.
##   wt_global = w_rails    * w
##   wt_sub    = w_subrails * w
########################################################################

wdata_clean <- km %>%
  filter(fh_wgs_ehr,
         !is.na(w), !is.na(w_rails), !is.na(w_subrails),
         !is.na(ehr_length), ehr_length > 0) %>%
  mutate(Age       = ifelse(Age > 89, 89, Age),
         wt_global = w_rails    * w,
         wt_sub    = w_subrails * w)

message("Analytic cohort (FH+WGS+EHR, both weights, ehr_length > 0): ", nrow(wdata_clean))
message("Breast-cancer events: ", sum(wdata_clean$has_breastcancerDX, na.rm = TRUE))

########################################################################
## Kaplan-Meier plots — one per weighting scheme.
## Strata: BRCA carrier status x any breast-cancer family history.
########################################################################

if (!dir.exists("../graph")) dir.create("../graph", recursive = TRUE)

surv_form <- Surv(start_age, end_age, has_breastcancerDX) ~ has_BRCA + has_BreastCancer_FH_any
km_theme  <- theme_bw() + theme(plot.title = element_text(hjust = 0.5))

## Short strata labels (order = survfit strata order: BRCA 0/1 x FH 0/1)
strata_labs <- c("No BRCA, No FH", "No BRCA, FH", "BRCA, No FH", "BRCA, FH")

build_km <- function(fit, title) {
  ggsurvplot(fit, data = wdata_clean, conf.int = TRUE, risk.table = FALSE,
             ggtheme = km_theme, legend = "bottom", legend.title = "",
             legend.labs = strata_labs, title = title,
             xlab = "Age", ylab = "Breast-cancer-free probability")
}

fit_unw    <- survfit(surv_form, data = wdata_clean)
fit_global <- survfit(surv_form, data = wdata_clean, weights = wt_global)
fit_sub    <- survfit(surv_form, data = wdata_clean, weights = wt_sub)

## survminer can't subset a formula that was passed as a variable (its stored
## call$formula is the symbol `surv_form`). Put the real formula back so
## ggsurvplot can read the strata.
fit_unw$call$formula    <- surv_form
fit_global$call$formula <- surv_form
fit_sub$call$formula    <- surv_form

## All three schemes in ONE figure, panels in a row, sharing a SINGLE legend
## at the bottom (ggarrange pulls the legend from the first panel).
km_panels <- list(build_km(fit_unw,    "Unweighted"),
                  build_km(fit_global, "Global RAILS"),
                  build_km(fit_sub,    "Subgroup RAILS"))
km_row <- ggpubr::ggarrange(plotlist = lapply(km_panels, `[[`, "plot"),
                            ncol = 3, nrow = 1,
                            common.legend = TRUE, legend = "bottom")
print(km_row)
ggsave("../graph/km_three_schemes_row.png", km_row, width = 16, height = 5.5, dpi = 300)

########################################################################
## Cumulative-risk (cumulative-incidence) plots — the complement of the KM
## curves. fun = "event" draws 1 - S(t), the cumulative probability of a
## breast-cancer diagnosis. Same strata, same three weighting schemes.
########################################################################

build_cuminc <- function(fit, title) {
  ggsurvplot(fit, data = wdata_clean, fun = "event",
             conf.int = TRUE, risk.table = FALSE,
             ggtheme = km_theme, legend = "bottom", legend.title = "",
             legend.labs = strata_labs, title = title,
             xlab = "Age", ylab = "Cumulative breast-cancer risk")
}

cuminc_panels <- list(build_cuminc(fit_unw,    "Unweighted"),
                      build_cuminc(fit_global, "Global RAILS"),
                      build_cuminc(fit_sub,    "Subgroup RAILS"))
cuminc_row <- ggpubr::ggarrange(plotlist = lapply(cuminc_panels, `[[`, "plot"),
                                ncol = 3, nrow = 1,
                                common.legend = TRUE, legend = "bottom")
print(cuminc_row)
ggsave("../graph/cuminc_three_schemes_row.png", cuminc_row, width = 16, height = 5.5, dpi = 300)

########################################################################
## Weight diagnostics — how much does each scheme actually reweight?
## If wt_global / wt_sub have tiny spread (CV) or w_rails ~ w_subrails
## (correlation ~ 1), then the curves SHOULD look alike — that's a data
## property, not a plotting error.
########################################################################

cat("\nWeight summaries (analytic cohort):\n")
print(round(sapply(wdata_clean[c("w", "w_rails", "w_subrails", "wt_global", "wt_sub")],
                   function(x) c(mean = mean(x), sd = sd(x), CV = sd(x) / mean(x),
                                 min = min(x), max = max(x))), 3))
cat("\ncor(w_rails, w_subrails) = ", round(cor(wdata_clean$w_rails, wdata_clean$w_subrails), 4),
    "\nmean |w_rails - w_subrails| / mean(w_rails) = ",
    round(mean(abs(wdata_clean$w_rails - wdata_clean$w_subrails)) / mean(wdata_clean$w_rails), 4), "\n")

########################################################################
## Overlaid comparison (pooled female cohort) so small gaps are visible
## side by side rather than across three separate panels.
########################################################################

ov_fit <- function(wt = NULL) {
  f <- if (is.null(wt)) survfit(Surv(start_age, end_age, has_breastcancerDX) ~ 1, data = wdata_clean)
       else            survfit(Surv(start_age, end_age, has_breastcancerDX) ~ 1, data = wdata_clean, weights = wt)
  tibble(time = f$time, surv = f$surv)
}

overlay <- bind_rows(
  ov_fit()                      %>% mutate(scheme = "Unweighted"),
  ov_fit(wdata_clean$wt_global) %>% mutate(scheme = "Global RAILS"),
  ov_fit(wdata_clean$wt_sub)    %>% mutate(scheme = "Subgroup RAILS (by sex)")
) %>%
  mutate(scheme = factor(scheme, levels = c("Unweighted", "Global RAILS", "Subgroup RAILS (by sex)")))

scheme_cols <- c("Unweighted" = "grey40",
                 "Global RAILS" = "#2166AC",
                 "Subgroup RAILS (by sex)" = "#B2182B")

p_ov <- ggplot(overlay, aes(time, surv, color = scheme)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = scheme_cols) +
  labs(x = "Age", y = "Breast-cancer-free probability", color = NULL,
       title = "Kaplan-Meier by weighting scheme (pooled female cohort)") +
  km_theme
print(p_ov)
ggsave("../graph/km_overlay_schemes.png", p_ov, width = 9, height = 6, dpi = 300)

## Cumulative-risk overlay: the same curves as 1 - S(t).
p_ov_risk <- ggplot(overlay, aes(time, 1 - surv, color = scheme)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = scheme_cols) +
  labs(x = "Age", y = "Cumulative breast-cancer risk", color = NULL,
       title = "Cumulative risk by weighting scheme (pooled female cohort)") +
  km_theme
print(p_ov_risk)
ggsave("../graph/cuminc_overlay_schemes.png", p_ov_risk, width = 9, height = 6, dpi = 300)

########################################################################
## (1) Weighted survival + 95% CI at ages 75/80/85/90 (pooled cohort).
## NOTE: with weights these are Greenwood CIs on weighted counts — treat as
## approximate (they can understate variance under heavy weighting).
########################################################################

surv_at <- function(wt, scheme, times = c(75, 80, 85, 90)) {
  f <- if (is.null(wt)) survfit(Surv(start_age, end_age, has_breastcancerDX) ~ 1, data = wdata_clean)
       else            survfit(Surv(start_age, end_age, has_breastcancerDX) ~ 1, data = wdata_clean, weights = wt)
  s <- summary(f, times = times, extend = TRUE)
  tibble(age = s$time, scheme = scheme,
         cell = sprintf("%.3f (%.3f, %.3f)", s$surv, s$lower, s$upper))
}

surv_table <- bind_rows(
  surv_at(NULL,                  "Unweighted"),
  surv_at(wdata_clean$wt_global, "Global RAILS"),
  surv_at(wdata_clean$wt_sub,    "Subgroup RAILS")
) %>%
  mutate(scheme = factor(scheme, levels = c("Unweighted", "Global RAILS", "Subgroup RAILS"))) %>%
  tidyr::pivot_wider(names_from = scheme, values_from = cell)

cat("\nBreast-cancer-free probability, surv (95% CI):\n")
print(as.data.frame(surv_table), row.names = FALSE)

########################################################################
## (2) Winsorize weights at the 99th percentile and re-draw the tail, so we
## can see whether the age>75 separation is real or driven by extreme weights.
########################################################################

cap <- function(x, p = 0.99) pmin(x, quantile(x, p, na.rm = TRUE))
wdata_clean <- wdata_clean %>%
  mutate(wt_global_tr = cap(wt_global), wt_sub_tr = cap(wt_sub))

tail_overlay <- bind_rows(
  ov_fit()                         %>% mutate(scheme = "Unweighted",     trim = "Untrimmed"),
  ov_fit(wdata_clean$wt_global)    %>% mutate(scheme = "Global RAILS",   trim = "Untrimmed"),
  ov_fit(wdata_clean$wt_sub)       %>% mutate(scheme = "Subgroup RAILS", trim = "Untrimmed"),
  ov_fit()                         %>% mutate(scheme = "Unweighted",     trim = "Winsorized (99th pct)"),
  ov_fit(wdata_clean$wt_global_tr) %>% mutate(scheme = "Global RAILS",   trim = "Winsorized (99th pct)"),
  ov_fit(wdata_clean$wt_sub_tr)    %>% mutate(scheme = "Subgroup RAILS", trim = "Winsorized (99th pct)")
) %>%
  mutate(scheme = factor(scheme, levels = c("Unweighted", "Global RAILS", "Subgroup RAILS")),
         trim   = factor(trim,   levels = c("Untrimmed", "Winsorized (99th pct)")))

p_tail <- ggplot(tail_overlay, aes(time, surv, color = scheme)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = c("Unweighted" = "grey40",
                                "Global RAILS" = "#2166AC",
                                "Subgroup RAILS" = "#B2182B")) +
  coord_cartesian(xlim = c(65, 102)) +
  facet_wrap(~trim) +
  labs(x = "Age", y = "Breast-cancer-free probability", color = NULL,
       title = "Tail (age > 65): untrimmed vs winsorized weights") +
  km_theme
print(p_tail)
ggsave("../graph/km_tail_trim_compare.png", p_tail, width = 12, height = 5.5, dpi = 300)

########################################################################
## Save the weighted analytic dataset
########################################################################

write_excel_csv(wdata_clean, "../data/weighted_data_brca_bysex.csv")
message("Wrote ../data/weighted_data_brca_bysex.csv")
