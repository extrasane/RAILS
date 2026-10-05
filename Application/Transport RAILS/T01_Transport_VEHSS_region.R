#### T01_Transport_VEHSS_region.R — Transportability of RAILS weights across
#### Census regions, validated against the CDC VEHSS regional benchmarks.
#### Run inside the AoU Researcher Workbench (Controlled Tier). PRELIMINARY run.
####
#### Notation follows the Subgroup RAILS manuscript; see DEFINITIONS.md.
#### For each target region a (target population U^a), three AoU SOURCE SAMPLES
#### are calibrated to the SAME targets T^{pop,a} (region a's ACS PUMS 40+ margins):
####   source  set        meaning                                   RAILS estimate
####   within  U_A^a      AoU participants living in a              S-RAILS (generalization)
####   all     U_A        all AoU participants                      P-RAILS (pooled)
####   out     U_A^{-a}   AoU participants living outside a         T-RAILS (transport)
#### and their weighted prevalences of AMD, glaucoma and cataract are compared
#### with the VEHSS estimate for region a (states aggregated by implied
#### population, exactly as in the VEHSS_comparison*.R scripts). The transport
#### contrast is delta_out_within = T-RAILS - S-RAILS (benchmark-free).
####
#### Nothing is queried here. The script reads the STAGE 1 / STAGE 2 outputs
#### that the VEHSS comparison scripts already wrote to ../data/:
####   vehss_aou_40plus.csv                 (VEHSS_comparison.R,          AMD)
####   vehss_glaucoma_aou_40plus.csv        (VEHSS_comparison_glaucoma.R, glaucoma)
####   vehss_cataract_aou_40plus.csv        (VEHSS_comparison_cataract.R, cataract)
####   vehss_pums_2022_bystate_40plus.csv   (stage 1 of any of them; identical 40+ PUMS)
#### plus the VEHSS_<AMD|GLAUCOMA|CATARACT>_{overall,combinations}.csv tables
#### built by VEHSS_data_eye.R. A phenotype whose AoU file is missing is skipped.
#### The covariates / cohort of the FIRST available AoU file define the weights;
#### the other phenotypes' case counts are joined to it by person_id (the three
#### stage-2 cohorts are the same 40+ query, so they should match exactly; the
#### overlap is reported).
####
#### Requires Sub_AoU_Fun.R and AoU_Fun.R (sourced exactly as in the
#### VEHSS_comparison*.R scripts: Sub_AoU_Fun.R with chdir = TRUE, which
#### sources the AoU_Fun.R sitting next to it). Locally, where the two files
#### live in different repo folders, AoU_Fun.R is sourced from
#### ../Global RAILS/RAILS Procedure/ first instead.
####
#### Estimator: fun.rails.threeway (two-way base + selected three-way terms,
#### LIFO raking) on the six covariates WITHOUT region (constant in the
#### target). ESTIMATOR = "auto" falls back to fun.rails.twoway for a source
#### whose three-way fit fails (typically a small within-region source); the
#### estimator used is recorded per (region, source). When no interaction is
#### selected, the RAILS weights are the raked base-model weights
#### (d_nps2_rake for three-way, d_nps1_rake for two-way), also recorded.
####
#### NOTE on the out source (and partly all): the pseudo-likelihood step assumes
#### the source is a subsample of the target population, which holds for within
#### only. For out the fitted "propensity" is a source-to-target density ratio
#### up to a constant; the raking step still enforces the target margins, so
#### the point estimate is well defined, but the analytic RAILS variance does
#### not apply. The CIs below are the naive linearized ones used throughout
#### the application scripts (05_Prevalence_Analysis.R / VEHSS_comparison.R);
#### they ignore calibration and selection, and are optimistic for every
#### source. Treat coverage as descriptive in this preliminary run.
####
#### Outputs (../data/transport_*.csv, ../graph/transport_*.png):
####   transport_fit_log.csv          estimator / terms / runtime per region x source
####   transport_overlap.csv          covariate overlap of each source with each target
####   transport_weight_diag.csv      ESS, weight variability, non-calibrated balance
####   transport_estimates.csv        prevalence vs VEHSS per region x source x method x outcome x age domain
####   transport_contrasts.csv        out - within and all - within per region x method x outcome x domain
####   transport_summary.csv          performance averaged over regions and outcomes
####   transport_heterogeneity.csv    regional outcome heterogeneity (marginal + conditional proxy)
####   transport_predictors.csv       one row per region x outcome x domain x source, all predictors of error
####   transport_size_matched.csv     (only if N_SIZE_MATCH > 0)
####   transport_cells.csv            cell-level weights (CONTAINS SMALL CELL COUNTS — keep in the Workbench)
####
#### AoU DISCLOSURE POLICY: estimates with fewer than SUPPRESS_MIN cases are
#### set to NA and their case count is masked. transport_cells.csv holds raw
#### cell counts and must not leave the Workbench.

suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(Matrix)
})

########################################################################
## Settings
########################################################################

SKIP_IF_CACHED <- FALSE   # TRUE = reuse transport_cells.csv / transport_fit_log.csv from an earlier run
DATA_DIR  <- "../data"
GRAPH_DIR <- "../graph"
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(GRAPH_DIR, showWarnings = FALSE, recursive = TRUE)

AGE_MIN        <- 40                    # VEHSS "40 years and older" (must match the stage-1/2 files)
AGE_BREAKS     <- c(40, 65, 85, Inf)    # VEHSS strata; re-cut here on both sides
MIN_CODE_DATES <- 2                     # distinct code dates required to be a case
ALPHA          <- 0.05                  # forward LRT threshold
ESTIMATOR      <- "auto"                # "threeway", "twoway", or "auto" (three-way, two-way on failure)
SUPPRESS_MIN   <- 20                    # AoU small-cell threshold (cases per estimate)
TRIM_Q         <- 0.99                  # RAILS weights capped at this weighted quantile (sensitivity method)
OVERLAP_ORDER  <- 2                     # membership model for the overlap c-statistic: 1 = main effects, 2 = + two-way
N_SIZE_MATCH   <- 0                     # > 0: draws of the all / out sources subsampled to the within-region size (re-runs RAILS each draw; slow)
SEED           <- 20261001

## Phenotypes: stage-2 AoU file, its case-date column, and the VEHSS tables
PHENOS <- tribble(
  ~pheno,     ~label,      ~aou_file,                        ~dates_col,          ~vehss_overall,                 ~vehss_combo,
  "amd",      "AMD",       "vehss_aou_40plus.csv",           "n_amd_dates",       "VEHSS_AMD_overall.csv",        "VEHSS_AMD_combinations.csv",
  "glaucoma", "Glaucoma",  "vehss_glaucoma_aou_40plus.csv",  "n_glaucoma_dates",  "VEHSS_GLAUCOMA_overall.csv",   "VEHSS_GLAUCOMA_combinations.csv",
  "cataract", "Cataract",  "vehss_cataract_aou_40plus.csv",  "n_cataract_dates",  "VEHSS_CATARACT_overall.csv",   "VEHSS_CATARACT_combinations.csv"
)
## 40+ PUMS by state (stage 1); the first file found is used
PUMS_FILES <- c("vehss_pums_2022_bystate_40plus.csv", "vehss_glaucoma_pums_2022_bystate_40plus.csv",
                "vehss_cataract_pums_2022_bystate_40plus.csv")

## Age domains: AoU agegroup levels and the VEHSS age rows that make up the benchmark
DOMAINS <- list(
  `40+`   = list(aou = c("40-64", "65-84", "85+"), vehss = "40+"),
  `40-64` = list(aou = "40-64",                    vehss = "40-64"),
  `65+`   = list(aou = c("65-84", "85+"),          vehss = c("65-84", "85+"))
)
## Benchmark comparability. Cataract VEHSS uses a Medicare-beneficiary
## denominator: 40-64 beneficiaries are disability / ESRD only, so only 65+
## is comparable with AoU's general population.
BENCH_NOT_COMPARABLE <- tribble(
  ~pheno,     ~domain,
  "cataract", "40+",
  "cataract", "40-64"
)

## Weight columns compared (all produced by fun.rails.threeway / twoway)
METHODS <- c(unweighted = "d_unweighted",   # equal weights: the source's raw prevalence
             rake1      = "d_cal1",         # one-way raking
             rake2      = "d_cal2",         # two-way raking
             nps2_rake  = "d_nps2_rake",    # two-way NPS raked to two-way margins
             rails      = "d_rails",        # RAILS (the method of interest)
             rails_trim = "d_rails_trim")   # RAILS capped at TRIM_Q, rescaled to the target total
SOURCES <- c("within", "all", "out")

########################################################################
## RAILS functions — same sourcing as the VEHSS_comparison*.R scripts
########################################################################

## Sub_AoU_Fun.R / AoU_Fun.R are looked up in the working directory, its
## parent, ../data, the repo folders, then anywhere below the parent / home.
find_fun_file <- function(name) {
  cands <- file.path(c(".", "..", DATA_DIR, "../Subgroup RAILS", "../Global RAILS/RAILS Procedure",
                       "~", "~/workspace", "/home/jupyter/workspace"), name)
  f <- cands[file.exists(cands)]
  if (length(f) == 0) {
    roots <- unique(c(normalizePath("..", mustWork = FALSE), path.expand("~")))
    roots <- roots[dir.exists(roots)]
    f <- unlist(lapply(roots, function(r)
      list.files(r, pattern = paste0("^", gsub(".", "[.]", name, fixed = TRUE), "$"),
                 recursive = TRUE, full.names = TRUE)))
  }
  if (length(f) == 0)
    stop(name, " not found. Upload it (Sub_AoU_Fun.R from Application/Subgroup RAILS/, AoU_Fun.R from ",
         "Application/Global RAILS/RAILS Procedure/) into ", normalizePath(getwd()), " and re-run.")
  normalizePath(f[1])
}
SUB_FUN_FILE <- find_fun_file("Sub_AoU_Fun.R")
if (!any(grepl("fun.rails.twoway", readLines(SUB_FUN_FILE, warn = FALSE), fixed = TRUE)))
  stop(SUB_FUN_FILE, " does not define fun.rails.twoway — upload the current Sub_AoU_Fun.R.")

if (file.exists(file.path(dirname(SUB_FUN_FILE), "AoU_Fun.R"))) {
  ## Workbench layout: both files side by side, Sub_AoU_Fun.R sources AoU_Fun.R
  source(SUB_FUN_FILE, chdir = TRUE)
  message("RAILS functions: ", SUB_FUN_FILE, " (+ AoU_Fun.R next to it)")
} else {
  ## Repo layout: AoU_Fun.R lives in Global RAILS/RAILS Procedure. Source it,
  ## then evaluate Sub_AoU_Fun.R without its own source("AoU_Fun.R") call.
  AOU_FUN_FILE <- find_fun_file("AoU_Fun.R")
  source(AOU_FUN_FILE)
  sub_exprs <- parse(SUB_FUN_FILE)
  is_source <- vapply(sub_exprs, function(e) is.call(e) && identical(e[[1]], as.name("source")), logical(1))
  for (e in sub_exprs[!is_source]) eval(e, envir = globalenv())
  message("RAILS functions: ", AOU_FUN_FILE, " + ", SUB_FUN_FILE)
}
stopifnot(exists("fun.rails.threeway"), exists("fun.rails.twoway"),
          exists("fun.nps"), exists("create_v3"))

########################################################################
## Shared definitions (identical to the VEHSS_comparison*.R scripts)
########################################################################

NAMES_6 <- c("agegroup", "sex", "edu", "homeown", "income", "race_eth")   # region is constant in a target

AGE_LABELS <- {
  b <- AGE_BREAKS
  k <- length(b) - 1
  c(paste0(b[seq_len(k - 1)], "-", b[seq_len(k - 1) + 1] - 1), paste0(b[k], "+"))
}
make_agegroup <- function(age) {
  as.character(cut(age, breaks = AGE_BREAKS, right = FALSE, labels = AGE_LABELS))
}
if (!identical(AGE_LABELS, c("40-64", "65-84", "85+"))) {
  warning("AGE_BREAKS are not the VEHSS strata: only the 40+ domain is compared.", immediate. = TRUE)
  DOMAINS <- list(`40+` = list(aou = AGE_LABELS, vehss = "40+"))
}

## State -> region map: official Census regions (as PUMS REGION and the VEHSS_comparison*.R scripts)
south     <- c("AL","AR","DC","DE","FL","GA","KY","LA","MD","MS","NC","OK","SC","TN","TX","VA","WV")
midwest   <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI")
northeast <- c("CT","ME","MA","NH","NJ","NY","PA","RI","VT")
west      <- c("AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY")
state_region_map <- c(
  setNames(rep("South",     length(south)),     south),
  setNames(rep("Midwest",   length(midwest)),   midwest),
  setNames(rep("Northeast", length(northeast)), northeast),
  setNames(rep("West",      length(west)),      west)
)
REGION_LEVELS <- c("Northeast", "Midwest", "South", "West")

harmonize_factors <- function(df) harmonize_levels(df) %>% na.omit()
harmonize_levels <- function(df) {
  df %>%
    mutate(
      agegroup = factor(as.character(agegroup), levels = AGE_LABELS),
      sex      = factor(sex,      levels = c("Female", "Male")),
      race_eth = factor(race_eth, levels = c("Hispanic", "NH Asian", "NH Black", "NH White", "Others")),
      income   = factor(income,   levels = c("<35k", "35k-50k", "50k-75k", "75k-100k", ">100k")),
      edu      = factor(edu,      levels = c("Less than highschool", "Some highschool",
                                             "Highschool graduate", "Some college",
                                             "College graduate or advanced")),
      homeown  = factor(ifelse(homeown == "Other", "Others", as.character(homeown)),
                        levels = c("Own", "Rent", "Others")),
      region   = factor(region,   levels = REGION_LEVELS)
    )
}

find_input <- function(name) {
  cands <- unique(c(name, file.path(DATA_DIR, name), file.path("../Subgroup RAILS", name), file.path("../Sub_VEHSS", name)))
  f <- cands[file.exists(cands)]
  if (length(f) == 0) return(NA_character_)
  f[1]
}

out_file <- function(x) file.path(DATA_DIR, paste0("transport_", x, ".csv"))
tic <- function() proc.time()[["elapsed"]]

########################################################################
## Small helpers: weighted prevalence, quantile, AUC, aggregation
########################################################################

## Weighted prevalence with a linearized 95% CI (as VEHSS_comparison.R)
wt.prev <- function(y, w) {
  ok <- !is.na(y) & !is.na(w); y <- y[ok]; w <- w[ok]
  if (length(y) == 0 || sum(w) <= 0) return(c(est = NA_real_, se = NA_real_))
  p <- sum(w * y) / sum(w)
  v <- sum(w^2 * (y - p)^2) / sum(w)^2
  c(est = p, se = sqrt(v))
}

## Quantile of x where each value carries frequency / weight f
wquantile <- function(x, f, q) {
  o <- order(x); x <- x[o]; cf <- cumsum(f[o]) / sum(f)
  x[which(cf >= q)[1]]
}

## P(score of class 1 > score of class 0), ties counted 1/2; w1 / w0 are the
## class-1 / class-0 weights carried by each row (rows may carry both)
wauc <- function(score, w1, w0) {
  d <- tibble(s = score, w1 = w1, w0 = w0) %>%
    group_by(s) %>% summarise(w1 = sum(w1), w0 = sum(w0), .groups = "drop") %>% arrange(s)
  below0 <- cumsum(d$w0) - d$w0
  sum(d$w1 * (below0 + 0.5 * d$w0)) / (sum(d$w1) * sum(d$w0))
}

agg_cells <- function(d, vars, wcol = "weight") {
  d %>% group_by(across(all_of(vars))) %>%
    summarise(weight = sum(.data[[wcol]]), .groups = "drop") %>%
    as.data.frame()
}

########################################################################
## STEP 1 — Load the stage-1 / stage-2 files of the VEHSS scripts
########################################################################

PHENOS <- PHENOS %>% mutate(aou_path = vapply(aou_file, find_input, character(1)))
missing_ph <- PHENOS$pheno[is.na(PHENOS$aou_path)]
if (length(missing_ph))
  warning("No stage-2 AoU file for: ", paste(missing_ph, collapse = ", "),
          " — run the matching VEHSS_comparison*.R stage 2 first. Skipped here.", immediate. = TRUE)
PHENOS <- PHENOS %>% filter(!is.na(aou_path))
if (nrow(PHENOS) == 0) stop("No stage-2 AoU file found in ", DATA_DIR, ".")

## Base cohort = the first available phenotype file
aou <- read_csv(PHENOS$aou_path[1], show_col_types = FALSE) %>%
  mutate(person_id = as.numeric(person_id),
         agegroup  = make_agegroup(age),
         region    = unname(state_region_map[state])) %>%
  filter(age >= AGE_MIN) %>%
  select(person_id, state, region, age, all_of(NAMES_6)) %>%
  harmonize_factors()
message("AoU ", AGE_MIN, "+ base cohort (", basename(PHENOS$aou_path[1]), "): ", nrow(aou), " participants")

## Case indicators y_<pheno>, joined by person_id
for (i in seq_len(nrow(PHENOS))) {
  ph <- PHENOS[i, ]
  yy <- read_csv(ph$aou_path, show_col_types = FALSE,
                 col_select = c("person_id", ph$dates_col)) %>%
    mutate(person_id = as.numeric(person_id),
           !!paste0("y_", ph$pheno) := as.integer(coalesce(as.numeric(.data[[ph$dates_col]]), 0) >= MIN_CODE_DATES)) %>%
    select(person_id, starts_with("y_")) %>%
    distinct(person_id, .keep_all = TRUE)
  aou <- aou %>% left_join(yy, by = "person_id")
  n_miss <- sum(is.na(aou[[paste0("y_", ph$pheno)]]))
  message(sprintf("  %-9s cases (>= %d dates): %d | base participants not in this file: %d",
                  ph$pheno, MIN_CODE_DATES, sum(aou[[paste0("y_", ph$pheno)]], na.rm = TRUE), n_miss))
  if (n_miss > 0.01 * nrow(aou))
    warning(ph$pheno, ": ", n_miss, " base participants missing from ", basename(ph$aou_path),
            " — its estimates use the overlap only.", immediate. = TRUE)
}
aou <- aou %>% mutate(weight = 1)

## PUMS 40+, rescaled to the 40+ population of the mapped states (as stage 3)
pums_path <- vapply(PUMS_FILES, find_input, character(1))
pums_path <- pums_path[!is.na(pums_path)][1]
if (is.na(pums_path)) stop("No 40+ PUMS file (", paste(PUMS_FILES, collapse = " / "), ") in ", DATA_DIR)
pums40 <- read_csv(pums_path, show_col_types = FALSE) %>%
  filter(age >= AGE_MIN) %>%
  mutate(region = unname(state_region_map[state]))
n_pop40 <- sum(pums40$PWGTP[!is.na(pums40$region)], na.rm = TRUE)
pums <- pums40 %>%
  mutate(agegroup = make_agegroup(age)) %>%
  select(state, region, PWGTP, age, all_of(NAMES_6)) %>%
  harmonize_factors() %>%
  mutate(weight = PWGTP / sum(PWGTP) * n_pop40)
message("PUMS ", AGE_MIN, "+ (", basename(pums_path), "): ", nrow(pums), " records, population ",
        format(round(n_pop40), big.mark = ","))

aou_agg  <- agg_cells(aou,  c(NAMES_6, "region"))
pums_agg <- agg_cells(pums, c(NAMES_6, "region"))

########################################################################
## STEP 2 — VEHSS regional benchmarks per outcome x age domain
########################################################################

read_vehss <- function(f) {
  readr::read_csv(f, show_col_types = FALSE,
                  col_types = readr::cols(age = "c", sex = "c", race = "c", state = "c",
                                          .default = readr::col_guess())) %>%
    transmute(age, sex, race, state, prev = as.numeric(prev),
              cases = as.numeric(cases), pop = as.numeric(pop))
}

bench <- bind_rows(lapply(seq_len(nrow(PHENOS)), function(i) {
  ph <- PHENOS[i, ]
  f_o <- find_input(ph$vehss_overall); f_c <- find_input(ph$vehss_combo)
  if (is.na(f_o)) { warning("VEHSS table ", ph$vehss_overall, " not found — ", ph$pheno, " has no benchmark."); return(NULL) }
  v <- bind_rows(read_vehss(f_o), if (!is.na(f_c)) read_vehss(f_c)) %>%
    filter(sex == "Both", race == "All", state != "US", !is.na(prev), !is.na(cases), !is.na(pop)) %>%
    mutate(region = unname(state_region_map[state])) %>%
    filter(!is.na(region))
  bind_rows(lapply(names(DOMAINS), function(dm) {
    ages <- DOMAINS[[dm]]$vehss
    st <- v %>% filter(age %in% ages) %>% distinct(state, age, .keep_all = TRUE) %>%
      group_by(state, region) %>%
      filter(n_distinct(age) == length(ages)) %>%          # a state enters only with every band of the domain
      summarise(cases = sum(cases), pop = sum(pop), .groups = "drop")
    if (nrow(st) == 0) return(NULL)
    by_r <- st %>% group_by(region) %>%
      summarise(bench = sum(cases) / sum(pop), bench_cases = sum(cases), bench_pop = sum(pop),
                bench_states = n(), .groups = "drop")
    ## Marginal regional heterogeneity: region vs the other three combined
    by_r %>% mutate(bench_others = (sum(bench_cases) - bench_cases) / (sum(bench_pop) - bench_pop),
                    het_marginal = bench - bench_others,
                    pheno = ph$pheno, domain = dm)
  }))
})) %>%
  left_join(BENCH_NOT_COMPARABLE %>% mutate(comparable = FALSE), by = c("pheno", "domain")) %>%
  mutate(comparable = coalesce(comparable, TRUE))
message("VEHSS regional benchmarks: ", nrow(bench), " region x outcome x domain rows")

########################################################################
## STEP 3 — RAILS: each source calibrated to each target region's margins
########################################################################

## Source cells for one (target region, source), collapsed over region
source_cells <- function(r, source, ids = NULL) {
  d <- if (is.null(ids)) aou else aou %>% filter(person_id %in% ids)
  d <- switch(source,
              within    = d %>% filter(region == r),
              all   = d,
              out = d %>% filter(region != r))
  agg_cells(d, NAMES_6)
}
target_cells <- function(r) agg_cells(pums %>% filter(region == r), NAMES_6)

## One-, two-, three-way target margins (a superset of what either estimator needs)
target_totals <- function(tgt) {
  twovars   <- combn(NAMES_6, 2, FUN = function(x) paste(x, collapse = ":"))
  threevars <- combn(NAMES_6, 3, FUN = function(x) paste(x, collapse = ":"))
  f   <- formula(paste0("~", paste(c(NAMES_6, twovars, threevars), collapse = "+")))
  mat <- sparse.model.matrix(f, data = tgt, keep.order = TRUE)
  setNames(as.numeric(Matrix::crossprod(mat, tgt$weight)), colnames(mat))
}

## Fit one source to one target. Returns the cell table with per-individual
## d_* columns (incl. d_rails_trim) plus the estimator bookkeeping.
fit_transport <- function(src, tgt, pop_totals) {
  nsiz <- sum(tgt$weight)
  try_fit <- function(est) {
    t0 <- tic()
    res <- tryCatch(
      suppressWarnings(
        if (est == "threeway")
          fun.rails.threeway(src, tgt, pop_totals, names_univar = NAMES_6, alpha = ALPHA, nsiz = nsiz)
        else
          fun.rails.twoway(src, tgt, pop_totals, names_univar = NAMES_6, alpha = ALPHA, nsiz = nsiz)),
      error = function(e) { message("    ", est, " failed: ", conditionMessage(e)); NULL })
    if (is.null(res)) return(NULL)
    base_col <- if (est == "threeway") "d_nps2_rake" else "d_nps1_rake"
    note <- "selected model raked"
    if (all(is.na(res$d_rails))) {
      res$d_rails <- res[[base_col]]
      note <- paste0("no interaction selected/converged: base model (", base_col, ")")
    }
    if (all(is.na(res$d_rails))) return(NULL)
    attr(res, "info") <- list(estimator = est, rails_note = note, seconds = tic() - t0)
    res
  }
  ests <- switch(ESTIMATOR, auto = c("threeway", "twoway"), threeway = "threeway", twoway = "twoway")
  res <- NULL
  for (est in ests) { res <- try_fit(est); if (!is.null(res)) break }
  if (is.null(res)) return(NULL)
  info <- attr(res, "info")
  ## Trimmed RAILS: cap the per-individual weight at its TRIM_Q quantile, rescale to the target total
  cap <- wquantile(res$d_rails, res$weight, TRIM_Q)
  res <- res %>% mutate(d_rails_trim = pmin(d_rails, cap),
                        d_rails_trim = d_rails_trim * nsiz / sum(d_rails_trim * weight))
  attr(res, "info") <- c(info, list(trim_cap = cap))
  res
}

F_CELLS <- out_file("cells"); F_LOG <- out_file("fit_log")
if (SKIP_IF_CACHED && file.exists(F_CELLS) && file.exists(F_LOG)) {
  message("SKIP_IF_CACHED = TRUE: reusing ", F_CELLS)
  ## levels only (no na.omit: benchmark weight columns may be NA); region borrowed from the target
  ## Files from runs before the renaming carry `cohort` = A_local / B_pooled / C_external
  old_to_new <- function(d) {
    if ("cohort" %in% names(d)) d <- d %>% rename(source = cohort)
    d %>% mutate(source = recode(source, A_local = "within", B_pooled = "all", C_external = "out"))
  }
  cells_all <- read_csv(F_CELLS, show_col_types = FALSE) %>% old_to_new() %>%
    mutate(region = target) %>% harmonize_levels() %>% select(-region)
  fit_log <- read_csv(F_LOG, show_col_types = FALSE) %>% old_to_new()
} else {
  cells_list <- list(); log_list <- list()
  for (r in REGION_LEVELS) {
    tgt <- target_cells(r); pt <- target_totals(tgt)
    for (co in SOURCES) {
      src <- source_cells(r, co)
      message(format(Sys.time(), "%H:%M:%S"), "  target ", r, " <- ", co, " (n = ", sum(src$weight), ")")
      res <- fit_transport(src, tgt, pt)
      if (is.null(res)) {
        warning("No weights for target ", r, " / ", co, ".", immediate. = TRUE)
        log_list[[length(log_list) + 1]] <- tibble(target = r, source = co, n_src = sum(src$weight),
                                                   estimator = NA_character_, rails_note = "FAILED")
        next
      }
      info <- attr(res, "info")
      cells_list[[length(cells_list) + 1]] <- res %>%
        select(all_of(NAMES_6), weight, all_of(unname(METHODS))) %>%
        mutate(target = r, source = co)
      log_list[[length(log_list) + 1]] <- tibble(
        target = r, source = co, n_src = sum(src$weight), n_cells_src = nrow(src),
        pop_target = sum(tgt$weight), estimator = info$estimator, rails_note = info$rails_note,
        selected_terms = res$selected_terms[1], calibrated_terms = res$calibrated_terms[1],
        trim_cap = info$trim_cap, seconds = round(info$seconds, 1))
    }
  }
  cells_all <- bind_rows(cells_list)
  fit_log   <- bind_rows(log_list)
  write_excel_csv(cells_all, F_CELLS)
  write_excel_csv(fit_log,   F_LOG)
}
print(fit_log %>% select(target, source, n_src, estimator, rails_note, seconds), n = Inf)

## Per-individual weights of one (target, source)
source_individuals <- function(r, co) {
  d <- switch(co,
              within    = aou %>% filter(region == r),
              all   = aou,
              out = aou %>% filter(region != r))
  w <- cells_all %>% filter(target == r, source == co) %>% select(all_of(NAMES_6), all_of(unname(METHODS)))
  d %>% left_join(w, by = NAMES_6)
}

########################################################################
## STEP 4 — Overlap (before weighting) and weight diagnostics
########################################################################

AGE5_BREAKS <- c(seq(AGE_MIN, 85, by = 5), Inf)
age5 <- function(a) cut(a, AGE5_BREAKS, right = FALSE)

overlap_list <- list(); diag_list <- list()
for (r in REGION_LEVELS) {
  tgt  <- target_cells(r)
  ptgt <- pums %>% filter(region == r)
  tgt_age_mean <- sum(ptgt$weight * ptgt$age) / sum(ptgt$weight)
  tgt_age_sd   <- sqrt(sum(ptgt$weight * (ptgt$age - tgt_age_mean)^2) / sum(ptgt$weight))
  tgt_age5     <- ptgt %>% group_by(band = age5(age)) %>% summarise(p_t = sum(weight), .groups = "drop") %>%
    mutate(p_t = p_t / sum(p_t))

  for (co in SOURCES) {
    src <- source_cells(r, co)

    ## (a) support: target population in joint cells without any source participant
    j <- tgt %>% rename(pop = weight) %>% left_join(src %>% rename(n = weight), by = NAMES_6)
    support_joint0 <- sum(j$pop[is.na(j$n)]) / sum(j$pop)
    ## (b) support: target share in two-way cells with < 5 source participants (mean / worst pair)
    pairs <- combn(NAMES_6, 2, simplify = FALSE)
    sh2 <- vapply(pairs, function(p) {
      t2 <- agg_cells(tgt, p) %>% rename(pop = weight)
      s2 <- agg_cells(src, p) %>% rename(n = weight)
      x  <- t2 %>% left_join(s2, by = p) %>% mutate(n = coalesce(n, 0))
      sum(x$pop[x$n < 5]) / sum(x$pop)
    }, numeric(1))
    ## (c) membership c-statistic: target (PUMS) vs source (AoU), class weights balanced
    mem <- bind_rows(tgt %>% mutate(z = 1, w = weight / sum(weight)),
                     src %>% mutate(z = 0, w = weight / sum(weight)))
    terms_m <- if (OVERLAP_ORDER >= 2)
      c(NAMES_6, combn(NAMES_6, 2, FUN = function(x) paste(x, collapse = ":"))) else NAMES_6
    auc <- tryCatch({
      fit <- suppressWarnings(glm(reformulate(terms_m, response = "z"), data = mem,
                                  family = quasibinomial(), weights = w))
      s <- predict(fit, mem, type = "link")
      wauc(s, w1 = mem$w * mem$z, w0 = mem$w * (1 - mem$z))
    }, error = function(e) NA_real_)

    overlap_list[[length(overlap_list) + 1]] <- tibble(
      target = r, source = co, n_src = sum(src$weight),
      support_joint0 = support_joint0, support_2way_lt5_mean = mean(sh2),
      support_2way_lt5_worst = max(sh2),
      support_2way_worst_pair = paste(pairs[[which.max(sh2)]], collapse = ":"),
      membership_auc = auc)

    ## Weight diagnostics per method
    ind <- source_individuals(r, co)
    wc  <- cells_all %>% filter(target == r, source == co)
    for (m in names(METHODS)) {
      w <- ind[[METHODS[m]]]
      if (is.null(w) || all(is.na(w))) next
      ok <- !is.na(w)
      sw <- sum(w[ok]); sw2 <- sum(w[ok]^2)
      ## Non-calibrated balance: detailed age (5-year bands, mean age) and the full 6-way joint
      a5 <- tibble(band = age5(ind$age[ok]), w = w[ok]) %>% group_by(band) %>%
        summarise(p_s = sum(w), .groups = "drop") %>% mutate(p_s = p_s / sum(p_s)) %>%
        full_join(tgt_age5, by = "band") %>% mutate(across(c(p_s, p_t), ~ coalesce(.x, 0)))
      src_age_mean <- sum(w[ok] * ind$age[ok]) / sw
      jc <- wc %>% transmute(across(all_of(NAMES_6)), tot = .data[[METHODS[m]]] * weight) %>%
        filter(!is.na(tot)) %>%
        full_join(tgt %>% rename(pop = weight), by = NAMES_6) %>%
        mutate(across(c(tot, pop), ~ coalesce(.x, 0)))
      diag_list[[length(diag_list) + 1]] <- tibble(
        target = r, source = co, method = m,
        n = sum(ok), n_weight_na = sum(!ok), sum_w = sw,
        ess = sw^2 / sw2, ess_pct = 100 * sw^2 / sw2 / sum(ok),
        cv_w = sd(w[ok]) / mean(w[ok]),
        max_over_median = max(w[ok]) / median(w[ok]),
        p99_over_p1 = unname(quantile(w[ok], 0.99) / quantile(w[ok], 0.01)),
        smd_age = (src_age_mean - tgt_age_mean) / tgt_age_sd,
        tvd_age5 = 0.5 * sum(abs(a5$p_s - a5$p_t)),
        tvd_joint6 = 0.5 * sum(abs(jc$tot / sum(jc$tot) - jc$pop / sum(jc$pop))))
    }
  }
}
overlap <- bind_rows(overlap_list)
wdiag   <- bind_rows(diag_list)
write_excel_csv(overlap, out_file("overlap"))
write_excel_csv(wdiag,   out_file("weight_diag"))
message("\n--- Overlap with target (before weighting) ---")
print(overlap %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf)

########################################################################
## STEP 5 — Weighted prevalence vs VEHSS
########################################################################

est_list <- list()
for (r in REGION_LEVELS) for (co in SOURCES) {
  ind <- source_individuals(r, co)
  for (ph in PHENOS$pheno) for (dm in names(DOMAINS)) {
    yv  <- ind[[paste0("y_", ph)]]
    sel <- ind$agegroup %in% DOMAINS[[dm]]$aou & !is.na(yv)
    for (m in names(METHODS)) {
      w <- ind[[METHODS[m]]][sel]
      if (all(is.na(w))) next
      pe <- wt.prev(yv[sel], w)
      est_list[[length(est_list) + 1]] <- tibble(
        target = r, source = co, method = m, pheno = ph, domain = dm,
        n = sum(sel & !is.na(ind[[METHODS[m]]])), n_cases = sum(yv[sel][!is.na(w)]),
        est = pe[["est"]], se = pe[["se"]])
    }
  }
}
estimates <- bind_rows(est_list) %>%
  mutate(suppress = n_cases < SUPPRESS_MIN,
         est = ifelse(suppress, NA_real_, est), se = ifelse(suppress, NA_real_, se),
         n_cases = ifelse(suppress, NA_integer_, as.integer(n_cases)),
         lb = est - qnorm(0.975) * se, ub = est + qnorm(0.975) * se) %>%
  left_join(bench %>% select(region, pheno, domain, bench, bench_states, comparable),
            by = c("target" = "region", "pheno", "domain")) %>%
  mutate(diff      = est - bench,                     # D
         abs_diff  = abs(diff),
         abs_rel_diff = abs_diff / bench,             # ARD
         pr        = est / bench,                     # prevalence ratio
         pr_lb     = lb / bench, pr_ub = ub / bench,  # benchmark has no aggregable SE: treated as fixed
         diff_lb   = lb - bench, diff_ub = ub - bench,
         covered   = !is.na(est) & !is.na(bench) & bench >= lb & bench <= ub,
         source    = factor(source, SOURCES), method = factor(method, names(METHODS)),
         target    = factor(target, REGION_LEVELS))
write_excel_csv(estimates, out_file("estimates"))

########################################################################
## STEP 6 — Regional outcome heterogeneity (does Y | X differ in vs out of r?)
## Unweighted cell-level binomial fits within AoU, per target x outcome x domain:
##   LRT of region-in/out x covariate interactions (conditional heterogeneity)
##   proxy = E_target[ m_out(X) - m_in(X) ]: the outcome-model part of the
##           transport bias of C, standardized to region r's PUMS margins
##   auc_y = discrimination of Y by the covariates (covariate-outcome strength)
########################################################################

het_list <- list()
for (r in REGION_LEVELS) for (ph in PHENOS$pheno) for (dm in names(DOMAINS)) {
  ycol <- paste0("y_", ph); ages <- DOMAINS[[dm]]$aou
  cl <- aou %>% filter(agegroup %in% ages, !is.na(.data[[ycol]])) %>%
    mutate(inr = as.integer(region == r)) %>%
    group_by(across(all_of(NAMES_6)), inr) %>%
    summarise(n = n(), cases = sum(.data[[ycol]]), .groups = "drop") %>%
    mutate(non = n - cases)
  vars <- NAMES_6[vapply(NAMES_6, function(v) n_distinct(cl[[v]]) > 1, logical(1))]
  out <- tryCatch(suppressWarnings({
    f0 <- reformulate(c(vars, "inr"), response = "cbind(cases, non)")
    f1 <- reformulate(c(vars, "inr", paste0(vars, ":inr")), response = "cbind(cases, non)")
    g0 <- glm(f0, data = cl, family = binomial()); g1 <- glm(f1, data = cl, family = binomial())
    fx <- reformulate(vars, response = "cbind(cases, non)")
    gA <- glm(fx, data = cl %>% filter(inr == 1), family = binomial())
    gC <- glm(fx, data = cl %>% filter(inr == 0), family = binomial())
    tg <- pums %>% filter(region == r, agegroup %in% ages) %>% agg_cells(NAMES_6)
    pA <- predict(gA, tg, type = "response"); pC <- predict(gC, tg, type = "response")
    tibble(lrt = g0$deviance - g1$deviance, lrt_df = g1$rank - g0$rank,
           lrt_p = pchisq(g0$deviance - g1$deviance, g1$rank - g0$rank, lower.tail = FALSE),
           std_in  = sum(tg$weight * pA) / sum(tg$weight),
           std_out = sum(tg$weight * pC) / sum(tg$weight),
           het_proxy = std_out - std_in,
           auc_y = wauc(predict(g0, cl, type = "link"), w1 = cl$cases, w0 = cl$non))
  }), error = function(e) tibble(lrt = NA_real_, het_proxy = NA_real_, auc_y = NA_real_))
  het_list[[length(het_list) + 1]] <- out %>% mutate(target = r, pheno = ph, domain = dm, .before = 1)
}
heterogeneity <- bind_rows(het_list) %>%
  left_join(bench %>% select(region, pheno, domain, bench, bench_others, het_marginal),
            by = c("target" = "region", "pheno", "domain"))
write_excel_csv(heterogeneity, out_file("heterogeneity"))

########################################################################
## STEP 7 — Synthesis: transport contrast (out - within), all - within, summaries, predictors
########################################################################

contrasts <- estimates %>%
  select(target, pheno, domain, method, source, est, abs_diff, comparable) %>%
  pivot_wider(names_from = source, values_from = c(est, abs_diff), names_expand = TRUE) %>%   # all 3 sources even if one failed
  mutate(delta_out_within   = est_out - est_within,      # benchmark-free
         delta_all_within   = est_all  - est_within,
         d_absdiff_out_within = abs_diff_out - abs_diff_within,   # > 0: T-RAILS further from VEHSS than S-RAILS
         d_absdiff_all_within = abs_diff_all  - abs_diff_within,
         target = as.character(target))
write_excel_csv(contrasts, out_file("contrasts"))

## Averages use only (target, outcome, domain, method) keys where ALL three
## sources have an un-suppressed estimate, so the sources are compared on the
## same set of comparisons (within, the smallest source, is suppressed most often)
summary_tab <- estimates %>%
  filter(comparable, !is.na(est), !is.na(bench)) %>%
  group_by(target, pheno, domain, method) %>%
  filter(n_distinct(source) == length(SOURCES)) %>%
  group_by(method, source) %>%
  summarise(n_comparisons    = n(),
            mean_abs_diff    = mean(abs_diff),
            mean_abs_rel_diff = mean(abs_rel_diff),
            median_pr        = median(pr),
            coverage_pct     = 100 * mean(covered),
            .groups = "drop") %>%
  left_join(contrasts %>% filter(comparable) %>% group_by(method) %>%
              summarise(out_worse_than_within_pct = 100 * mean(d_absdiff_out_within > 0, na.rm = TRUE),
                        mean_d_absdiff_out_within  = mean(d_absdiff_out_within, na.rm = TRUE),
                        all_worse_than_within_pct = 100 * mean(d_absdiff_all_within > 0, na.rm = TRUE),
                        mean_d_absdiff_all_within  = mean(d_absdiff_all_within, na.rm = TRUE), .groups = "drop"),
            by = "method")
write_excel_csv(summary_tab, out_file("summary"))
message("\n--- Performance vs VEHSS, averaged over regions x outcomes x domains (comparable benchmarks) ---")
print(summary_tab %>% mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf, width = Inf)

## One row per (target, outcome, domain, source): the RAILS error and its candidate predictors
predictors <- estimates %>%
  filter(method == "rails") %>%
  select(target, source, pheno, domain, comparable, est, bench, diff, abs_diff, abs_rel_diff) %>%
  left_join(estimates %>% filter(method == "unweighted") %>%
              select(target, source, pheno, domain, est_unweighted = est),
            by = c("target", "source", "pheno", "domain")) %>%
  mutate(calibration_shift = abs(est - est_unweighted),
         target = as.character(target), source = as.character(source)) %>%
  left_join(overlap %>% select(target, source, membership_auc, support_joint0, support_2way_lt5_mean),
            by = c("target", "source")) %>%
  left_join(wdiag %>% filter(method == "rails") %>% select(target, source, ess, ess_pct, cv_w, tvd_joint6),
            by = c("target", "source")) %>%
  left_join(heterogeneity %>% select(target, pheno, domain, het_marginal, het_proxy, lrt_p, auc_y),
            by = c("target", "pheno", "domain"))
write_excel_csv(predictors, out_file("predictors"))

## Exploratory association of the RAILS absolute error with each predictor
## (Spearman, comparable benchmarks; only ~4 regions x 3 outcomes x 3 domains)
pred_cols <- c("membership_auc", "support_joint0", "ess", "cv_w", "tvd_joint6", "bench",
               "het_marginal", "het_proxy", "auc_y", "calibration_shift")
assoc <- predictors %>% filter(comparable, !is.na(abs_diff)) %>%
  summarise(across(all_of(pred_cols), ~ suppressWarnings(cor(abs(.x), abs_diff, method = "spearman",
                                                             use = "complete.obs")))) %>%
  pivot_longer(everything(), names_to = "predictor", values_to = "spearman_with_abs_error")
message("\n--- Spearman correlation of |RAILS - VEHSS| with |predictor| (exploratory) ---")
print(assoc %>% mutate(spearman_with_abs_error = round(spearman_with_abs_error, 2)))

########################################################################
## STEP 8 (optional) — size-matched all / out sources: subsample to the within-region size
########################################################################

if (N_SIZE_MATCH > 0) {
  set.seed(SEED)
  sm_list <- list()
  for (r in REGION_LEVELS) {
    tgt <- target_cells(r); pt <- target_totals(tgt)
    n_A <- sum(aou$region == r)
    for (co in c("all", "out")) {
      pool <- if (co == "all") aou$person_id else aou$person_id[aou$region != r]
      for (b in seq_len(N_SIZE_MATCH)) {
        ids <- sample(pool, n_A)
        res <- fit_transport(source_cells(r, co, ids), tgt, pt)
        if (is.null(res)) next
        ind <- aou %>% filter(person_id %in% ids) %>%
          left_join(res %>% select(all_of(NAMES_6), d_unweighted, d_rails), by = NAMES_6)
        for (ph in PHENOS$pheno) for (dm in names(DOMAINS)) {
          yv <- ind[[paste0("y_", ph)]]; sel <- ind$agegroup %in% DOMAINS[[dm]]$aou & !is.na(yv)
          for (m in c("unweighted", "rails")) {
            pe <- wt.prev(yv[sel], ind[[METHODS[m]]][sel])
            sm_list[[length(sm_list) + 1]] <- tibble(target = r, source = co, draw = b, method = m,
                                                     pheno = ph, domain = dm, est = pe[["est"]],
                                                     n_cases = sum(yv[sel]))
          }
        }
      }
    }
  }
  size_matched <- bind_rows(sm_list) %>%
    group_by(target, source, method, pheno, domain) %>%
    summarise(draws = n(), est_mean = mean(est), est_sd = sd(est),
              min_cases = min(n_cases), .groups = "drop") %>%
    mutate(est_mean = ifelse(min_cases < SUPPRESS_MIN, NA_real_, est_mean),
           est_sd   = ifelse(min_cases < SUPPRESS_MIN, NA_real_, est_sd)) %>%
    select(-min_cases) %>%
    left_join(bench %>% select(region, pheno, domain, bench, comparable),
              by = c("target" = "region", "pheno", "domain")) %>%
    mutate(abs_diff_mean = abs(est_mean - bench))
  write_excel_csv(size_matched, out_file("size_matched"))
}

########################################################################
## STEP 9 — Figures (../graph/transport_fig*.png)
## One encoding throughout:
##   colour = source sample (within blue, all orange, out aqua — validated
##            categorical slots 1-3), filled dot = RAILS, open dot = unweighted,
##   black tick = VEHSS benchmark, rows = target region,
##   panels = outcome x age domain (comparable benchmarks only).
## Signed differences use a blue-grey-red diverging scale centred at 0.
##   fig1  estimates vs VEHSS, unweighted -> RAILS per source
##   fig2  ratio to VEHSS on one shared log scale
##   fig3  heatmap of the benchmark-free contrasts T-RAILS - S-RAILS and P-RAILS - S-RAILS
##   fig4  method x source mean error (dumbbell)
##   fig5  overlap and weight diagnostics per region x source
##   fig6  transport contrast (T-RAILS - S-RAILS) vs the concept-shift proxy
##   fig7  RAILS error vs ESS and vs overlap
##   fig8  size-matched out-of-region source (only if N_SIZE_MATCH > 0)
########################################################################

SOURCE_COLS   <- c(within = "#2a78d6", all = "#eb6834", out = "#1baf7a")
## Legend labels: set notation of the manuscript (plotmath). The colour marks the
## SOURCE SAMPLE for every method; S-/P-/T-RAILS name only its RAILS estimate.
SOURCE_LABELS <- c(within = expression("Within-region " * (U[A]^a)),
                   all    = expression("All-region " * (U[A])),
                   out    = expression("Out-of-region " * (U[A]^{-a})))
CAPTION_EST <- "RAILS on the within-, all- and out-of-region source = S-RAILS, P-RAILS, T-RAILS."
SOURCE_OFFSET <- c(within = 0.22, all = 0, out = -0.22)   # within on top within a row
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; SURFACE <- "#fcfcfb"
DIV_LOW <- "#2a78d6"; DIV_MID <- "#f0efec"; DIV_HIGH <- "#e34948"
METHOD_LABELS <- c(unweighted = "Unweighted", rake1 = "One-way raking", rake2 = "Two-way raking",
                   nps2_rake = "Two-way NPS + raking", rails = "RAILS", rails_trim = "RAILS, trimmed")

pheno_lab    <- setNames(PHENOS$label, PHENOS$pheno)
panel_levels <- unlist(lapply(PHENOS$pheno, function(p) paste0(pheno_lab[p], " ", names(DOMAINS))))
mk_panel <- function(pheno, domain) factor(paste0(pheno_lab[pheno], " ", domain), levels = panel_levels)
row_y    <- function(target) match(as.character(target), rev(REGION_LEVELS))   # Northeast on top
pp       <- function(x) 100 * x

theme_transport <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(plot.background  = element_rect(fill = SURFACE, colour = NA),
          panel.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.text  = element_text(colour = INK2), axis.title = element_text(colour = INK2),
          strip.text = element_text(colour = INK, face = "bold", hjust = 0),
          plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = INK2), plot.caption = element_text(colour = INK2, hjust = 0),
          legend.position = "bottom", legend.title = element_text(colour = INK2),
          legend.text = element_text(colour = INK))
}
scale_source <- function() scale_colour_manual(values = SOURCE_COLS, breaks = SOURCES,
                                               labels = SOURCE_LABELS[SOURCES],
                                               name = "Source sample", drop = FALSE)
scale_region_y <- function() scale_y_continuous(breaks = seq_along(REGION_LEVELS), labels = rev(REGION_LEVELS),
                                                expand = expansion(add = 0.5))
save_fig <- function(g, name, w, h) {
  f <- file.path(GRAPH_DIR, paste0("transport_", name, ".png"))
  ggsave(f, g, width = w, height = h, dpi = 200, bg = SURFACE)
  message("Wrote ", f)
}
CAPTION_CI <- paste0(CAPTION_EST, "\nBars: naive linearized 95% CI (ignores calibration and model selection). Estimates with < 20 cases are suppressed.")

## Comparable estimates with plotting coordinates
est_plot <- estimates %>%
  filter(comparable) %>%
  mutate(panel = mk_panel(pheno, domain),
         y     = row_y(target) + SOURCE_OFFSET[as.character(source)])
bench_plot <- est_plot %>% distinct(panel, target, bench) %>% mutate(y = row_y(target))

## ---- fig1: estimates vs VEHSS, the unweighted -> RAILS move per source ----
f1_r <- est_plot %>% filter(method == "rails")
f1_u <- est_plot %>% filter(method == "unweighted")
f1_m <- f1_u %>% select(panel, target, source, y, x0 = est) %>%
  inner_join(f1_r %>% select(panel, target, source, x1 = est), by = c("panel", "target", "source"))
g1 <- ggplot() +
  geom_segment(data = bench_plot, aes(x = pp(bench), xend = pp(bench), y = y - 0.4, yend = y + 0.4),
               colour = INK, linewidth = 0.9) +
  geom_segment(data = f1_m, aes(x = pp(x0), xend = pp(x1), y = y, yend = y, colour = source),
               linewidth = 0.4, alpha = 0.6) +
  geom_point(data = f1_u, aes(x = pp(est), y = y, colour = source), shape = 21, fill = SURFACE,
             size = 2, stroke = 0.8) +
  geom_linerange(data = f1_r, aes(xmin = pp(lb), xmax = pp(ub), y = y, colour = source), linewidth = 0.6) +
  geom_point(data = f1_r, aes(x = pp(est), y = y, colour = source), size = 2.3) +
  scale_source() + scale_region_y() +
  facet_wrap(~ panel, scales = "free_x", ncol = 3, drop = TRUE) +
  labs(x = "Prevalence (%)", y = NULL,
       title = "Prevalence by source sample, all calibrated to the target region's margins",
       subtitle = "Black tick: VEHSS. Open dot: unweighted source sample; filled dot + bar: RAILS. The line shows how far calibration moved the estimate.",
       caption = CAPTION_CI) +
  theme_transport()
save_fig(g1, "fig1_estimates", 11, 2 + 2.2 * ceiling(nlevels(droplevels(est_plot$panel)) / 3))

## ---- fig2: ratio to VEHSS on a shared log scale (all outcomes comparable) ----
f2 <- est_plot %>% filter(method %in% c("unweighted", "rails"), !is.na(pr)) %>%
  mutate(prow = as.numeric(factor(panel, levels = rev(levels(droplevels(panel))))),
         yy   = prow + SOURCE_OFFSET[as.character(source)],
         pr_lb = pmax(pr_lb, 1 / 512), pr_ub = pmax(pr_ub, 1 / 512))   # a CI reaching <= 0 is cut at the axis floor
prow_labels <- rev(levels(droplevels(f2$panel)))
g2 <- ggplot(f2, aes(y = yy, colour = source)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.5) +
  geom_linerange(data = f2 %>% filter(method == "rails"), aes(xmin = pr_lb, xmax = pr_ub), linewidth = 0.6) +
  geom_point(data = f2 %>% filter(method == "unweighted"), aes(x = pr), shape = 21, fill = SURFACE,
             size = 1.8, stroke = 0.7) +
  geom_point(data = f2 %>% filter(method == "rails"), aes(x = pr), size = 2.2) +
  scale_source() +
  scale_x_continuous(trans = "log2", breaks = 4^(-4:1),
                     labels = c("1/256", "1/64", "1/16", "1/4", "1", "4")) +
  scale_y_continuous(breaks = seq_along(prow_labels), labels = prow_labels, expand = expansion(add = 0.5)) +
  facet_wrap(~ factor(target, REGION_LEVELS), nrow = 1) +
  labs(x = "Prevalence ratio, AoU / VEHSS (log scale; 1 = agreement)", y = NULL,
       title = "Agreement with VEHSS by target region",
       subtitle = "Filled dot + bar: RAILS; open dot: unweighted. Left of 1 = AoU lower than VEHSS.",
       caption = CAPTION_CI) +
  theme_transport()
save_fig(g2, "fig2_ratio_to_vehss", 11, 1.5 + 0.55 * length(prow_labels))

## ---- fig3 / fig3b: benchmark-free contrasts T-RAILS - S-RAILS and P-RAILS - S-RAILS ----
## fig3  in percentage points; fig3b relative to S-RAILS (%), so a small
## contrast on a rare outcome (AMD) is not washed out by cataract's scale
f3 <- contrasts %>% filter(method == "rails", comparable) %>%
  transmute(target, pheno, domain,
            abs_out = pp(delta_out_within), abs_all = pp(delta_all_within),
            rel_out = 100 * delta_out_within / est_within, rel_all = 100 * delta_all_within / est_within) %>%
  pivot_longer(c(abs_out, abs_all, rel_out, rel_all), names_to = c("scale", "contrast"),
               names_sep = "_", values_to = "shift") %>%
  mutate(contrast = factor(contrast, c("out", "all"),
                           c("T-RAILS - S-RAILS (transport contrast)",
                             "P-RAILS - S-RAILS")),
         panel  = mk_panel(pheno, domain),
         target = factor(target, REGION_LEVELS))
heat_fig <- function(sc, fmt, legend_name, subtitle, file) {
  d   <- f3 %>% filter(scale == sc) %>% mutate(lab = ifelse(is.na(shift), "–", sprintf(fmt, shift)))
  lim <- max(abs(d$shift), 1e-6, na.rm = TRUE)
  g <- ggplot(d, aes(x = target, y = fct_rev(droplevels(panel)), fill = shift)) +
    geom_tile(colour = SURFACE, linewidth = 1) +
    geom_text(aes(label = lab, colour = abs(shift) > 0.6 * lim), size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "#ffffff"), na.value = INK2) +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH, midpoint = 0,
                         limits = c(-lim, lim), na.value = GRID, name = legend_name) +
    facet_wrap(~ contrast, nrow = 1) +
    labs(x = "Target region", y = NULL,
         title = "How much does the source sample change the RAILS estimate?",
         subtitle = subtitle) +
    theme_transport() + theme(panel.grid.major = element_blank(), legend.position = "right")
  save_fig(g, file, 10, 1.5 + 0.45 * nlevels(droplevels(d$panel)))
}
heat_fig("abs", "%+.2f", "Shift (pp)",
         "Percentage points; no benchmark involved. Near 0 = out-of-region participants reproduce S-RAILS. – = suppressed.",
         "fig3_transport_shift_heatmap")
heat_fig("rel", "%+.0f%%", "Shift (% of S-RAILS)",
         "Relative to S-RAILS, so outcomes of different prevalence are comparable. – = suppressed.",
         "fig3b_transport_shift_relative")

## ---- fig4: method x source mean relative error, per outcome (dumbbell) ----
## Relative error per outcome: pooled over outcomes, AMD's large benchmark gap
## (EHR-diagnosed vs modeled) would dominate every average.
f4_keys <- estimates %>%
  filter(comparable, !is.na(est), !is.na(bench)) %>%
  group_by(target, pheno, domain, method) %>%
  filter(n_distinct(source) == length(SOURCES)) %>% ungroup()
f4 <- bind_rows(f4_keys %>% mutate(outcome = pheno_lab[pheno]),
                f4_keys %>% mutate(outcome = "All outcomes")) %>%
  group_by(outcome, method, source) %>%
  summarise(value = 100 * mean(abs_rel_diff), n = n(), .groups = "drop") %>%
  mutate(outcome = factor(outcome, c(PHENOS$label, "All outcomes")),
         method  = factor(METHOD_LABELS[as.character(method)], rev(METHOD_LABELS)))
f4_rng <- f4 %>% group_by(outcome, method) %>% summarise(lo = min(value), hi = max(value), .groups = "drop")
f4_missing <- setdiff(METHOD_LABELS, as.character(unique(f4$method)))
g4 <- ggplot(f4, aes(y = method)) +
  geom_segment(data = f4_rng, aes(x = lo, xend = hi, yend = method), colour = GRID, linewidth = 2) +
  geom_point(aes(x = value, colour = source), size = 2.8) +
  scale_source() +
  facet_wrap(~ outcome, scales = "free_x", nrow = 1) +
  labs(x = "Mean |AoU - VEHSS| / VEHSS (%)", y = NULL,
       title = "Error against VEHSS by weighting method and source sample",
       subtitle = "Averaged over regions x age domains where all three source samples have an estimate (same comparisons for each).",
       caption = paste0(CAPTION_EST, if (length(f4_missing))
         paste0("\nNot shown (no comparison with all three sources converged/un-suppressed): ",
                paste(f4_missing, collapse = ", "), ". See transport_weight_diag.csv (n_weight_na).") else "")) +
  theme_transport()
save_fig(g4, "fig4_method_summary", 12, 4.5)

## ---- fig5: overlap and weight diagnostics per region x source ----
d_rails <- wdiag %>% filter(method == "rails")
d_unw   <- wdiag %>% filter(method == "unweighted")
f5 <- bind_rows(
  overlap %>% transmute(target, source, measure = "Target-vs-source c-statistic (AUC), before weighting\n0.5 = same covariate mix, 1 = fully separable",
                        before = NA_real_, after = membership_auc),
  overlap %>% transmute(target, source, measure = "% of target population (PUMS) in X-cells\nwith NO participant in the AoU source",
                        before = NA_real_, after = 100 * support_joint0),
  d_rails %>% transmute(target, source, measure = "Effective sample size, RAILS\n(log10)",
                        before = NA_real_, after = log10(ess)),
  d_rails %>% transmute(target, source, measure = "Weight CV, RAILS", before = NA_real_, after = cv_w),
  d_unw %>% select(target, source, before = tvd_joint6) %>%
    inner_join(d_rails %>% select(target, source, after = tvd_joint6), by = c("target", "source")) %>%
    mutate(measure = "Joint 6-way cell mismatch (TVD)\nunweighted -> RAILS"),
  d_unw %>% transmute(target, source, before = abs(smd_age)) %>%
    inner_join(d_rails %>% transmute(target, source, after = abs(smd_age)), by = c("target", "source")) %>%
    mutate(measure = "|Mean-age SMD|\nunweighted -> RAILS")
) %>%
  mutate(measure = factor(measure, unique(measure)),
         y = row_y(target) + SOURCE_OFFSET[as.character(source)],
         source = factor(source, SOURCES))
f5_ref <- tibble(measure = factor("|Mean-age SMD|\nunweighted -> RAILS", levels(f5$measure)), x = 0.1)
g5 <- ggplot(f5, aes(y = y, colour = source)) +
  geom_vline(data = f5_ref, aes(xintercept = x), colour = INK2, linewidth = 0.4) +   # usual 0.1 imbalance cut
  geom_segment(data = f5 %>% filter(!is.na(before)), aes(x = before, xend = after, yend = y),
               linewidth = 0.4, alpha = 0.6, show.legend = FALSE,
               arrow = arrow(length = unit(0.08, "inches"), type = "closed")) +
  geom_point(data = f5 %>% filter(!is.na(before)), aes(x = before), shape = 21, fill = SURFACE,
             size = 1.8, stroke = 0.7) +
  geom_point(aes(x = after), size = 2.3) +
  scale_source() + scale_region_y() +
  facet_wrap(~ measure, scales = "free_x", ncol = 3) +
  labs(x = NULL, y = NULL,
       title = "Overlap with the target and what the weighting cost",
       subtitle = "Balance panels (bottom row): open dot = unweighted, arrow tip = RAILS; these features are NOT calibrated, so lower after weighting is a real gain.",
       caption = paste0(CAPTION_EST,
         "\nc-statistic: AUC of a logistic model (X main effects + two-way terms) separating the target region's PUMS records",
         "\n    (weights d_i^B) from the AoU source participants, each class weighted to total 1; computed before any weighting.",
         "\nX-cell: one of the 2,250 combinations of the 6 covariates. Empty = the combination occurs in the region's PUMS",
         "\n    (target population > 0) but has no participant in the AoU source sample, so calibration cannot represent it directly.",
         "\nGrey line in the age panel: SMD = 0.1, the usual imbalance threshold.")) +
  theme_transport()
save_fig(g5, "fig5_diagnostics", 11, 7.2)

## ---- fig6: transport contrast (T-RAILS - S-RAILS) vs the concept-shift proxy ----
f6 <- contrasts %>% filter(method == "rails", comparable) %>%
  select(target, pheno, domain, delta_out_within) %>%
  inner_join(heterogeneity %>% select(target, pheno, domain, het_proxy, lrt_p),
             by = c("target", "pheno", "domain")) %>%
  filter(!is.na(delta_out_within), !is.na(het_proxy)) %>%
  mutate(sig = ifelse(!is.na(lrt_p) & lrt_p < 0.05, "p < 0.05", "p >= 0.05"),
         region_lab = c(Northeast = "NE", Midwest = "MW", South = "S", West = "W")[target],
         outcome = factor(pheno_lab[pheno], PHENOS$label),
         domain  = factor(domain, names(DOMAINS)))
## Per-outcome symmetric square limits (free scales, diagonal stays at 45 degrees)
f6_lim <- f6 %>% group_by(outcome) %>%
  summarise(l = max(abs(pp(c(delta_out_within, het_proxy))), 0.05, na.rm = TRUE) * 1.15, .groups = "drop") %>%
  reframe(x = c(-l, l), y = c(-l, l), .by = outcome)
label_fun <- if (requireNamespace("ggrepel", quietly = TRUE)) ggrepel::geom_text_repel else geom_text
g6 <- ggplot(f6, aes(x = pp(het_proxy), y = pp(delta_out_within))) +
  geom_blank(data = f6_lim, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, colour = GRID, linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = GRID, linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4) +
  geom_point(aes(shape = domain, fill = sig), colour = INK, size = 2.6, stroke = 0.7) +
  label_fun(aes(label = region_lab), size = 2.8, colour = INK2) +
  scale_shape_manual(values = c(`40+` = 21, `40-64` = 22, `65+` = 24), name = "Age domain") +
  scale_fill_manual(values = c(`p < 0.05` = INK, `p >= 0.05` = SURFACE),
                    limits = c("p < 0.05", "p >= 0.05"), drop = FALSE,
                    name = "In/out covariate effects differ (LRT)") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = INK, size = 2.6, stroke = 0.7))) +
  facet_wrap(~ outcome, scales = "free") +
  labs(x = "Concept-shift proxy: mean over target of [m_out(X) - m_within(X)] (pp)",
       y = "Transport contrast: T-RAILS - S-RAILS (pp)",
       title = "Is the transport contrast explained by concept shift (outcome given X differs across regions)?",
       subtitle = "On the diagonal: the contrast is what the regional difference in outcome given X predicts. Off it near the y-axis: look at balance and overlap (fig5).",
       caption = "Axes differ by outcome (each panel symmetric around 0). Labels: target region.") +
  theme_transport() + theme(aspect.ratio = 1)
save_fig(g6, "fig6_shift_vs_heterogeneity", 11, 5)

## ---- fig7: RAILS error vs ESS and vs overlap, one row per outcome ----
## Rows by outcome with free y: pooled, the outcome-level benchmark gap forms
## horizontal bands that hide any within-outcome trend.
f7 <- predictors %>% filter(comparable, !is.na(abs_diff)) %>%
  transmute(source = factor(source, SOURCES), outcome = factor(pheno_lab[pheno], PHENOS$label),
            domain = factor(domain, names(DOMAINS)), rel_err = 100 * abs_rel_diff,
            `Effective sample size (log10)` = log10(ess),
            `Target-vs-source c-statistic (AUC; 0.5 = same covariate mix)` = membership_auc) %>%
  pivot_longer(-c(source, outcome, domain, rel_err), names_to = "predictor", values_to = "value")
g7 <- ggplot(f7, aes(x = value, y = rel_err, colour = source, shape = domain)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_source() +
  scale_shape_manual(values = c(`40+` = 16, `40-64` = 15, `65+` = 17), name = "Age domain") +
  facet_grid(outcome ~ predictor, scales = "free") +
  labs(x = NULL, y = "|RAILS - VEHSS| / VEHSS (%)",
       title = "Does the error grow as precision or overlap falls?",
       subtitle = "One dot per region x age domain x source; rows = outcome (separate y scales).",
       caption = CAPTION_EST) +
  theme_transport() + theme(strip.text.y = element_text(angle = 0))
save_fig(g7, "fig7_error_vs_ess_overlap", 10, 2 + 2.2 * nrow(PHENOS))

## ---- fig8 (optional): S-RAILS vs T-RAILS (full and size-matched out-of-region source) ----
if (N_SIZE_MATCH > 0 && exists("size_matched")) {
  V8 <- c("S-RAILS (within-region)", "T-RAILS (out-of-region, full)", "T-RAILS (out-of-region, size-matched)")
  f8 <- bind_rows(
    est_plot %>% filter(method == "rails", source %in% c("within", "out")) %>%
      transmute(panel, target = as.character(target), version = ifelse(source == "within", V8[1], V8[2]),
                est, sd = NA_real_),
    size_matched %>% filter(method == "rails", source == "out", comparable) %>%
      transmute(panel = mk_panel(pheno, domain), target, version = V8[3],
                est = est_mean, sd = est_sd)) %>%
    mutate(version = factor(version, V8),
           y = row_y(target) + c(0.22, 0, -0.22)[as.integer(version)])
  g8 <- ggplot(f8, aes(y = y)) +
    geom_segment(data = bench_plot, aes(x = pp(bench), xend = pp(bench), y = y - 0.4, yend = y + 0.4),
                 colour = INK, linewidth = 0.9) +
    geom_linerange(aes(xmin = pp(est - sd), xmax = pp(est + sd), colour = version), linewidth = 0.6) +
    geom_point(aes(x = pp(est), colour = version, shape = version), size = 2.3, fill = SURFACE, stroke = 0.8) +
    scale_colour_manual(values = setNames(SOURCE_COLS[c("within", "out", "out")], V8), name = NULL) +
    scale_shape_manual(values = setNames(c(16, 16, 21), V8), name = NULL) +
    scale_region_y() +
    facet_wrap(~ panel, scales = "free_x", ncol = 3, drop = TRUE) +
    labs(x = "RAILS prevalence (%)", y = NULL,
         title = "Is the out-of-region source's advantage only its sample size?",
         subtitle = sprintf("Open dot + bar: mean +/- SD over %d draws of the out-of-region source subsampled to the within-region size. Black tick: VEHSS.", N_SIZE_MATCH)) +
    theme_transport()
  save_fig(g8, "fig8_size_matched", 11, 2 + 2.2 * ceiling(nlevels(droplevels(f8$panel)) / 3))
}

message("\nDone. Tables: ", DATA_DIR, "/transport_*.csv | figures: ", GRAPH_DIR, "/transport_fig*.png",
        "\n(transport_cells.csv contains raw cell counts — keep it inside the Workbench.)")
