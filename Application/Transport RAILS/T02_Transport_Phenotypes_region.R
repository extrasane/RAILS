#### T02_Transport_Phenotypes_region.R — Transportability of RAILS weights for
#### the GBD phenotype panel (prevalence AND association), by Census region.
#### Run inside the AoU Researcher Workbench (the SUBGROUP workspace), AFTER
#### 09_Sub_RAILS_region.R and Build_Disease_Indicators.R (and, for the
#### sex-specific exclusion list, 10_State_Prevalence_Maps.R with
#### EXCLUDE_SEX_SPECIFIC <- TRUE). Notation: DEFINITIONS.md (manuscript).
####
#### There is NO external benchmark here: every comparison is between source
#### samples calibrated to the SAME target, so the results describe how much
#### the choice of source changes the estimates, not which one is correct.
####
#### For each target region a (target population U^a, targets T^{pop,a} from
#### the manuscript's PUMS cell file dt_agg_pums_v2.csv, adults 18+, original
#### S-RAILS age groups 18-24 / 25-44 / 45-64 / 65-74 / 75+), three AoU source
#### samples are calibrated with RAILS to T^{pop,a}:
####   source  set        meaning                                 RAILS estimate
####   within  U_A^a      AoU participants living in a            S-RAILS
####   all     U_A        all AoU participants                    P-RAILS
####   out     U_A^{-a}   AoU participants living outside a       T-RAILS
#### G-RAILS (the manuscript's global weights w_rails, averaged over U_A^a) is
#### carried as a reference.
####
#### STEP 1  Calibration: 4 regions x 3 sources, fun.rails.threeway (two-way
####         fallback), six covariates (region constant in the target).
#### STEP 2  Check: the within-region weights must REPRODUCE w_subrails of
####         09_Sub_RAILS_region.R (same cells, targets and function).
#### PART A  Prevalence of every phenotype j in region a:
####         unweighted / one-way raking / RAILS for each source, plus G-RAILS;
####         contrasts T-RAILS - S-RAILS (transport contrast) and P-RAILS -
####         S-RAILS (pp, relative, logit); Jensen-Shannon divergence d^a and
####         the manuscript's D_j, E_j and leverage region for S-T, S-P and
####         S-G (unit = region, L = 4); sign consistency of the transport
####         contrast across phenotypes within each region; concept-shift proxy
####         and LRT of region-in/out x covariate interactions.
#### PART B  Association of every phenotype with the demographic covariates:
####         weighted logistic regression, CONDITIONAL (phenotype ~ all six
####         covariates) and MARGINAL (phenotype ~ one of sex / age group /
####         race-ethnicity), per source and weighting; transport contrast in
####         beta (T - S, P - S) with an approximate z; comparison with the
####         prevalence contrast on the same log-odds scale.
####
#### Computation. Every weight is constant within a covariate cell, so all
#### estimates (prevalences, linearized SEs, logistic MLEs, sandwich SEs) are
#### computed EXACTLY from cell-level case counts — no person-level loops.
#### Standard errors are the naive linearized / sandwich ones (they ignore
#### calibration and model selection), as in T01.
####
#### Outputs (../data/transport2_*.csv, ../graph/transport2_fig*.png):
####   transport2_fit_log.csv            estimator / terms per region x source
####   transport2_srails_check.csv       within-region weights vs w_subrails (09)
####   transport2_weight_diag.csv        n, ESS, weight CV per region x source
####   transport2_prevalence.csv         region x source x method x phenotype
####   transport2_contrasts.csv          S, P, T, G per region x phenotype + contrasts, JS
####   transport2_DE.csv                 D_j, E_j, leverage region for S-T, S-P, S-G
####   transport2_sign_consistency.csv   share of phenotypes with T-RAILS > S-RAILS per region
####   transport2_concept_shift.csv      concept-shift proxy + LRT per region x phenotype
####   transport2_beta.csv               beta, SE per region x source x method x model x phenotype x term
####   transport2_beta_contrasts.csv     beta(T) - beta(S), beta(P) - beta(S), z, + prevalence logit contrast
####   transport2_beta_summary.csv       summary of |contrasts| by model x term x region
####   transport2_cells.csv              cell-level weights (CONTAINS SMALL CELL COUNTS — keep in the Workbench)
####
#### AoU DISCLOSURE POLICY: prevalences with fewer than SUPPRESS_MIN cases are
#### NA; a beta is NA when its exposure level or the reference level has fewer
#### than SUPPRESS_MIN cases in the source. transport2_cells.csv must not leave
#### the Workbench.

suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(Matrix)
})

########################################################################
## Settings
########################################################################

SKIP_IF_CACHED    <- FALSE   # TRUE = reuse transport2_cells.csv / transport2_fit_log.csv
FETCH_FROM_BUCKET <- FALSE   # TRUE = gsutil cp the three pipeline files from $WORKSPACE_BUCKET/data/ (as 09 does)
DATA_DIR  <- "../data"
GRAPH_DIR <- "../graph"
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(GRAPH_DIR, showWarnings = FALSE, recursive = TRUE)

## Inputs (looked up in the working directory first, then DATA_DIR)
F_PUMS_AGG <- "dt_agg_pums_v2.csv"                          # PUMS cells, 7 covariates (manuscript targets)
F_AOU_AGG  <- "dt_agg_aou_v2.csv"                           # AoU cells, 7 covariates
F_AOU_RAW  <- "aou_raking_dt.csv"                           # AoU individuals (person_id + covariates)
F_DISEASE  <- "raking_wts_w_diseases_2_code_requirement.csv" # person_id, state, region, w_rails, phenotypes
F_SUBW     <- "dt_sub_aou_region.csv"                       # person_id, w_subrails (09), for the check
F_SEXSPEC  <- "sex_specific_phenotypes_excluded_nosex.csv"  # Disease column (10, _nosex run)

EXCLUDE_SEX_SPECIFIC <- TRUE   # drop sex-specific phenotypes (needed for the sex association)
ALPHA          <- 0.05         # forward LRT threshold (as 09)
ESTIMATOR      <- "auto"       # "threeway", "twoway", or "auto" (three-way, two-way on failure)
SUPPRESS_MIN   <- 20           # AoU small-cell threshold
MIN_CASES_BETA <- 100          # Part B: a phenotype needs >= this many cases in the source
RUN_CONCEPT_SHIFT <- TRUE      # Part A concept-shift proxy + LRT (~4 GLMs per region x phenotype)
RUN_BETA          <- TRUE      # Part B
BETA_MAX_ABS      <- 15        # |beta| above this (separation) is set to NA

## Reference levels for the association models
REF_LEVELS <- c(agegroup = "45-64", sex = "Female", race_eth = "NH White")

NAMES_6 <- c("agegroup", "sex", "edu", "homeown", "income", "race_eth")   # as 09_Sub_RAILS_region.R
REGION_LEVELS <- c("Northeast", "Midwest", "South", "West")
SOURCES <- c("within", "all", "out")
METHODS <- c(unweighted = "d_unweighted", rake1 = "d_cal1", rails = "d_rails")

########################################################################
## RAILS functions — same sourcing as T01 / the VEHSS_comparison*.R scripts
########################################################################

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
  source(SUB_FUN_FILE, chdir = TRUE)                       # Workbench layout
  message("RAILS functions: ", SUB_FUN_FILE, " (+ AoU_Fun.R next to it)")
} else {
  AOU_FUN_FILE <- find_fun_file("AoU_Fun.R")               # repo layout
  source(AOU_FUN_FILE)
  sub_exprs <- parse(SUB_FUN_FILE)
  is_source <- vapply(sub_exprs, function(e) is.call(e) && identical(e[[1]], as.name("source")), logical(1))
  for (e in sub_exprs[!is_source]) eval(e, envir = globalenv())
  message("RAILS functions: ", AOU_FUN_FILE, " + ", SUB_FUN_FILE)
}
stopifnot(exists("fun.rails.threeway"), exists("fun.rails.twoway"), exists("fun.nps"), exists("create_v3"))

########################################################################
## Shared definitions — IDENTICAL to 09_Sub_RAILS_region.R
########################################################################

harmonize_levels <- function(df) {
  df %>%
    mutate(
      sex      = factor(sex,      levels = c("Female", "Male")),
      race_eth = factor(race_eth, levels = c("Hispanic", "NH Asian", "NH Black", "NH White", "Others")),
      income   = factor(income,   levels = c("<35k", "35k-50k", "50k-75k", "75k-100k", ">100k")),
      agegroup = factor(agegroup, levels = c("18-24", "25-44", "45-64", "65-74", "75+")),
      edu      = factor(edu,      levels = c("Less than highschool", "Some highschool",
                                             "Highschool graduate", "Some college",
                                             "College graduate or advanced")),
      homeown  = factor(ifelse(homeown == "Other", "Others", as.character(homeown)),
                        levels = c("Own", "Rent", "Others")),
      region   = factor(region,   levels = REGION_LEVELS)
    )
}
harmonize_factors <- function(df) harmonize_levels(df) %>% na.omit()

## State -> region map: official Census (PUMS REGION), identical to 09
south     <- c("AL","AR","DC","DE","FL","GA","KY","LA","MD","MS","NC","OK","SC","TN","TX","VA","WV")  # official Census (PUMS REGION): + DC, DE, MD, OK
midwest   <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI")
northeast <- c("CT","ME","MA","NH","NJ","NY","PA","RI","VT")  # official Census: DE, MD moved to South
west      <- c("AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY")
state_region_map <- c(
  setNames(rep("South",     length(south)),     south),
  setNames(rep("Midwest",   length(midwest)),   midwest),
  setNames(rep("Northeast", length(northeast)), northeast),
  setNames(rep("West",      length(west)),      west)
)

find_input <- function(name, required = TRUE) {
  cands <- unique(c(name, file.path(DATA_DIR, name)))
  f <- cands[file.exists(cands)]
  if (length(f) == 0) {
    if (required) stop(name, " not found in ", paste(cands, collapse = " / "))
    return(NA_character_)
  }
  f[1]
}
out_file <- function(x) file.path(DATA_DIR, paste0("transport2_", x, ".csv"))
tic <- function() proc.time()[["elapsed"]]

agg_cells <- function(d, vars, wcol = "weight") {
  d %>% group_by(across(all_of(vars))) %>%
    summarise(weight = sum(.data[[wcol]]), .groups = "drop") %>% as.data.frame()
}
wquantile <- function(x, f, q) { o <- order(x); x <- x[o]; cf <- cumsum(f[o]) / sum(f); x[which(cf >= q)[1]] }
H_bern  <- function(p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); -p * log(p) - (1 - p) * log(1 - p) }
js_bern <- function(p, q) H_bern((p + q) / 2) - 0.5 * H_bern(p) - 0.5 * H_bern(q)
logit   <- function(p) log(p / (1 - p))

########################################################################
## STEP 0 — Load the pipeline files
########################################################################

if (FETCH_FROM_BUCKET) {
  my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
  for (fname in c(F_PUMS_AGG, F_AOU_AGG, F_AOU_RAW))
    system(paste0("gsutil cp ", my_bucket, "/data/", fname, " ."), intern = TRUE)
}

pums_agg <- read_csv(find_input(F_PUMS_AGG), show_col_types = FALSE) %>% harmonize_factors()
aou_agg  <- read_csv(find_input(F_AOU_AGG),  show_col_types = FALSE) %>% harmonize_factors()

## AoU individuals, built exactly as in 09_Sub_RAILS_region.R
dt_aou <- read_csv(find_input(F_AOU_RAW), show_col_types = FALSE) %>%
  select(!any_of(c("gender", "careplace"))) %>%
  na.omit() %>%
  mutate(
    race_eth = case_when(
      ethnicity == "Yes"                       ~ "Hispanic",
      race      == "Asian"                     ~ "NH Asian",
      race      == "Black or African American" ~ "NH Black",
      race      == "White"                     ~ "NH White",
      TRUE                                     ~ "Others"),
    agegroup = case_when(
      age <= 24            ~ "18-24",
      age > 24 & age <= 44 ~ "25-44",
      age > 44 & age <= 64 ~ "45-64",
      age > 64 & age <= 74 ~ "65-74",
      TRUE                 ~ "75+"),
    region = unname(state_region_map[state])) %>%
  harmonize_factors() %>%
  mutate(person_id = as.numeric(person_id)) %>%
  select(person_id, region, all_of(NAMES_6))

## Consistency of the individual file with the AoU cell file, by region
chk <- dt_aou %>% count(region, name = "n_individuals") %>%
  left_join(aou_agg %>% group_by(region) %>% summarise(n_cells_file = sum(weight), .groups = "drop"),
            by = "region")
message("AoU individuals vs ", F_AOU_AGG, " (by region):"); print(chk)
if (any(abs(chk$n_individuals - chk$n_cells_file) > 0.01 * chk$n_cells_file, na.rm = TRUE))
  warning("The individual file and the AoU cell file differ by > 1% in some region: calibration uses ",
          F_AOU_AGG, ", estimation uses the individuals.", immediate. = TRUE)

## Phenotypes
dz <- read_csv(find_input(F_DISEASE), show_col_types = FALSE) %>% mutate(person_id = as.numeric(person_id))
disease_cols <- setdiff(names(dz), c("person_id", "state", "region", "w_rails", "w_subrails"))
if (EXCLUDE_SEX_SPECIFIC) {
  f_ss <- find_input(F_SEXSPEC, required = FALSE)
  if (is.na(f_ss)) {
    warning("EXCLUDE_SEX_SPECIFIC = TRUE but ", F_SEXSPEC, " not found — run 10_State_Prevalence_Maps.R ",
            "with EXCLUDE_SEX_SPECIFIC <- TRUE first. ALL phenotypes are kept; the sex association of ",
            "sex-specific phenotypes is then meaningless.", immediate. = TRUE)
  } else {
    ss <- read_csv(f_ss, show_col_types = FALSE)$Disease
    disease_cols <- setdiff(disease_cols, ss)
    message("Sex-specific phenotypes excluded: ", length(intersect(ss, names(dz))))
  }
}
message("Phenotypes analysed: ", length(disease_cols))

ind <- dt_aou %>%
  inner_join(dz %>% select(person_id, w_rails, all_of(disease_cols)), by = "person_id") %>%
  filter(!is.na(w_rails))
message("Participants with covariates, region, phenotypes and w_rails: ", nrow(ind),
        " (of ", nrow(dt_aou), " individuals)")

########################################################################
## STEP 1 — RAILS: each source sample calibrated to each region's targets
########################################################################

target_cells <- function(a) agg_cells(pums_agg %>% filter(region == a), NAMES_6)
source_cells <- function(a, src) {
  d <- switch(src,
              within = aou_agg %>% filter(region == a),
              all    = aou_agg,
              out    = aou_agg %>% filter(region != a))
  agg_cells(d, NAMES_6)
}
## One-, two-, three-way targets (as fun.sub.rails.threeway builds them)
target_totals <- function(tgt) {
  twovars   <- combn(NAMES_6, 2, FUN = function(x) paste(x, collapse = ":"))
  threevars <- combn(NAMES_6, 3, FUN = function(x) paste(x, collapse = ":"))
  f   <- formula(paste0("~", paste(c(NAMES_6, twovars, threevars), collapse = "+")))
  mat <- sparse.model.matrix(f, data = tgt, keep.order = TRUE)
  setNames(as.numeric(Matrix::crossprod(mat, tgt$weight)), colnames(mat))
}
fit_transport <- function(src, tgt, pop_totals) {
  nsiz <- sum(tgt$weight)
  try_fit <- function(est) {
    t0 <- tic()
    res <- tryCatch(suppressWarnings(
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
  res
}

F_CELLS <- out_file("cells"); F_LOG <- out_file("fit_log")
if (SKIP_IF_CACHED && file.exists(F_CELLS) && file.exists(F_LOG)) {
  message("SKIP_IF_CACHED = TRUE: reusing ", F_CELLS)
  cells_all <- read_csv(F_CELLS, show_col_types = FALSE) %>%
    mutate(region = target) %>% harmonize_levels() %>% select(-region)
  fit_log <- read_csv(F_LOG, show_col_types = FALSE)
} else {
  cells_list <- list(); log_list <- list()
  for (a in REGION_LEVELS) {
    tgt <- target_cells(a); pt <- target_totals(tgt)
    for (s in SOURCES) {
      src <- source_cells(a, s)
      message(format(Sys.time(), "%H:%M:%S"), "  target ", a, " <- ", s, " (n = ", sum(src$weight), ")")
      res <- fit_transport(src, tgt, pt)
      if (is.null(res)) {
        warning("No weights for target ", a, " / ", s, ".", immediate. = TRUE)
        log_list[[length(log_list) + 1]] <- tibble(target = a, source = s, n_src = sum(src$weight),
                                                   estimator = NA_character_, rails_note = "FAILED")
        next
      }
      info <- attr(res, "info")
      cells_list[[length(cells_list) + 1]] <- res %>%
        select(all_of(NAMES_6), weight, all_of(unname(METHODS))) %>%
        mutate(target = a, source = s)
      log_list[[length(log_list) + 1]] <- tibble(
        target = a, source = s, n_src = sum(src$weight), pop_target = sum(tgt$weight),
        estimator = info$estimator, rails_note = info$rails_note,
        selected_terms = res$selected_terms[1], calibrated_terms = res$calibrated_terms[1],
        seconds = round(info$seconds, 1))
    }
  }
  cells_all <- bind_rows(cells_list); fit_log <- bind_rows(log_list)
  write_excel_csv(cells_all, F_CELLS); write_excel_csv(fit_log, F_LOG)
}
print(fit_log %>% select(target, source, n_src, estimator, rails_note), n = Inf)

########################################################################
## Cell machinery. U6 = all 3,750 combinations of the six covariates
## (5 x 2 x 5 x 3 x 5 x 5). Person i -> (cell c, region r).
## n_r[c, r]: participants; Y_r[[r]][c, j]: cases of phenotype j.
########################################################################

lev <- lapply(setNames(NAMES_6, NAMES_6), function(v) levels(ind[[v]]))
U6  <- expand.grid(lev, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
key <- function(d) do.call(paste, c(lapply(NAMES_6, function(v) as.integer(factor(d[[v]], levels = lev[[v]]))), sep = "_"))
U6_key <- key(U6); NC <- nrow(U6)

ind_c <- match(key(ind), U6_key)
ind_r <- as.integer(ind$region)
Dmat  <- as(data.matrix(ind[disease_cols]), "dgCMatrix")      # persons x phenotypes
n_r   <- matrix(0, NC, length(REGION_LEVELS)); Y_r <- vector("list", length(REGION_LEVELS))
wG_r  <- matrix(NA_real_, NC, length(REGION_LEVELS))           # G-RAILS weight per (cell, region)
for (r in seq_along(REGION_LEVELS)) {
  idx <- which(ind_r == r)
  M   <- sparseMatrix(i = seq_along(idx), j = ind_c[idx], x = 1, dims = c(length(idx), NC))
  n_r[, r]  <- as.numeric(Matrix::colSums(M))
  Y_r[[r]]  <- as.matrix(Matrix::crossprod(M, Dmat[idx, , drop = FALSE]))   # cells x phenotypes
  wsum      <- as.numeric(Matrix::crossprod(M, ind$w_rails[idx]))
  wG_r[, r] <- ifelse(n_r[, r] > 0, wsum / n_r[, r], NA_real_)
}
colnames(n_r) <- REGION_LEVELS
rm(Dmat); gc()

## G-RAILS weights should be constant within (cell, region): report the spread
wg_sd <- tapply(ind$w_rails, paste(ind_c, ind_r), function(x) if (length(x) > 1) sd(x) / mean(x) else 0)
message("G-RAILS weight: max within-cell coefficient of variation = ", signif(max(wg_sd, na.rm = TRUE), 3),
        " (0 = constant within cell, as expected)")

## PUMS target population per cell and region; rho_a
pop6 <- sapply(REGION_LEVELS, function(a) {
  t <- target_cells(a); v <- numeric(NC); v[match(key(t), U6_key)] <- t$weight; v })
rho <- colSums(pop6) / sum(pop6)

member_regions <- function(a, s) switch(s, within = a, all = REGION_LEVELS, out = setdiff(REGION_LEVELS, a))
src_counts <- function(a, s) {
  R <- match(member_regions(a, s), REGION_LEVELS)
  list(n = rowSums(n_r[, R, drop = FALSE]), Y = Reduce(`+`, Y_r[R]))
}
## Per-cell weight vector (per person) for (target, source, method); 0 where the cell has no member
cell_weights <- function(a, s, m) {
  if (m == "grails") { w <- wG_r[, match(a, REGION_LEVELS)]; return(ifelse(is.na(w), 0, w)) }
  cw <- cells_all %>% filter(target == a, source == s)
  w <- numeric(NC)
  if (nrow(cw) == 0) return(rep(NA_real_, NC))
  w[match(key(cw), U6_key)] <- cw[[METHODS[m]]]
  w[is.na(w)] <- 0
  w
}

########################################################################
## STEP 2 — Check: within-region RAILS must reproduce w_subrails (09)
########################################################################

f_sub <- find_input(F_SUBW, required = FALSE)
if (!is.na(f_sub)) {
  subw <- read_csv(f_sub, show_col_types = FALSE) %>% mutate(person_id = as.numeric(person_id)) %>%
    select(person_id, w_subrails) %>% inner_join(dt_aou, by = "person_id") %>% filter(!is.na(w_subrails))
  sub_c <- match(key(subw), U6_key); sub_r <- as.integer(subw$region)
  srails_check <- bind_rows(lapply(seq_along(REGION_LEVELS), function(r) {
    a <- REGION_LEVELS[r]; k <- sub_r == r
    ours <- cell_weights(a, "within", "rails")[sub_c[k]]
    ref  <- subw$w_subrails[k]
    tibble(target = a, n = sum(k),
           max_abs_rel_diff = max(abs(ours - ref) / ref, na.rm = TRUE),
           median_abs_rel_diff = median(abs(ours - ref) / ref, na.rm = TRUE),
           cor = suppressWarnings(cor(ours, ref, use = "complete.obs")))
  }))
  write_excel_csv(srails_check, out_file("srails_check"))
  message("\n--- Check: within-region RAILS vs w_subrails of 09_Sub_RAILS_region.R ---")
  print(srails_check)
  if (any(srails_check$max_abs_rel_diff > 1e-4, na.rm = TRUE))
    warning("Within-region RAILS does NOT reproduce w_subrails (max relative difference > 1e-4). ",
            "Check that the same input files / Sub_AoU_Fun.R version were used before interpreting.",
            immediate. = TRUE)
} else {
  message(F_SUBW, " not found — S-RAILS reproduction check skipped.")
}

## Weight diagnostics (per person weights, from cells)
wdiag <- bind_rows(lapply(REGION_LEVELS, function(a) bind_rows(lapply(SOURCES, function(s) {
  cnt <- src_counts(a, s); w <- cell_weights(a, s, "rails"); ok <- cnt$n > 0 & w > 0
  sw <- sum(w[ok] * cnt$n[ok]); sw2 <- sum(w[ok]^2 * cnt$n[ok]); n <- sum(cnt$n[ok])
  mw <- sw / n; vw <- sum(cnt$n[ok] * (w[ok] - mw)^2) / (n - 1)
  tibble(target = a, source = s, n = n, ess = sw^2 / sw2, ess_pct = 100 * sw^2 / sw2 / n, cv_w = sqrt(vw) / mw,
         n_in_unweighted_cells = sum(cnt$n[cnt$n > 0 & !(w > 0)]))
}))))
write_excel_csv(wdiag, out_file("weight_diag"))
print(wdiag %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf)

########################################################################
## PART A — Prevalence of every phenotype under every source x method
########################################################################

prev_one <- function(a, s, m) {
  cnt <- src_counts(a, s); w <- cell_weights(a, s, m)
  if (m == "unweighted") w <- as.numeric(cnt$n > 0)
  if (all(is.na(w))) return(NULL)
  ok <- cnt$n > 0 & w > 0
  W  <- sum(w[ok] * cnt$n[ok]); W2 <- sum(w[ok]^2 * cnt$n[ok])
  num <- as.numeric(crossprod(cnt$Y[ok, , drop = FALSE], w[ok]))
  A   <- as.numeric(crossprod(cnt$Y[ok, , drop = FALSE], w[ok]^2))
  p   <- num / W
  v   <- ((1 - p)^2 * A + p^2 * (W2 - A)) / W^2          # sum w^2 (y - p)^2 / (sum w)^2
  tibble(target = a, source = s, method = m, pheno = disease_cols,
         n = sum(cnt$n[ok]), n_cases = as.numeric(colSums(cnt$Y[ok, , drop = FALSE])),
         est = p, se = sqrt(v))
}
prevalence <- bind_rows(c(
  lapply(REGION_LEVELS, function(a) bind_rows(lapply(SOURCES, function(s)
    bind_rows(lapply(names(METHODS), function(m) prev_one(a, s, m)))))),
  lapply(REGION_LEVELS, function(a) prev_one(a, "within", "grails")))) %>%
  mutate(suppress = n_cases < SUPPRESS_MIN,
         est = ifelse(suppress, NA_real_, est), se = ifelse(suppress, NA_real_, se),
         n_cases = ifelse(suppress, NA_real_, n_cases))
write_excel_csv(prevalence, out_file("prevalence"))

## Contrasts per region x phenotype (RAILS): S = within, P = all, T = out, G = G-RAILS
contrasts <- prevalence %>%
  filter((method == "rails") | (method == "grails")) %>%
  mutate(est_name = case_when(method == "grails" ~ "G", source == "within" ~ "S",
                              source == "all" ~ "P", source == "out" ~ "T")) %>%
  select(target, pheno, est_name, est) %>%
  pivot_wider(names_from = est_name, values_from = est) %>%
  left_join(prevalence %>% filter(method == "unweighted") %>% select(target, pheno, source, est) %>%
              pivot_wider(names_from = source, values_from = est, names_prefix = "unw_"),
            by = c("target", "pheno")) %>%
  mutate(delta_out_within  = T - S,                       # transport contrast (pp / 100)
         delta_all_within  = P - S,
         delta_grails_within = G - S,
         rel_out_within    = (T - S) / S,
         rel_all_within    = (P - S) / S,
         log2_ratio_TS     = log2(T / S),
         log2_ratio_PS     = log2(P / S),
         logit_delta_TS    = logit(T) - logit(S),
         js_ST = js_bern(S, T), js_SP = js_bern(S, P), js_SG = js_bern(S, G),
         rho = rho[target])
write_excel_csv(contrasts, out_file("contrasts"))

## D_j, E_j and leverage region (manuscript application metrics, unit = region, L = 4).
## rho renormalized over the regions that enter (un-suppressed), as in the VEHSS scripts.
de_one <- function(d, rho_a, regions) {
  ok <- is.finite(d) & is.finite(rho_a)
  if (sum(ok) < 2) return(tibble(n_regions = sum(ok), D = NA_real_, E = NA_real_,
                                 leverage_region = NA_character_, leverage_share = NA_real_))
  r  <- rho_a[ok] / sum(rho_a[ok]); c_a <- r * d[ok]; D <- sum(c_a)
  if (D <= 0) return(tibble(n_regions = sum(ok), D = D, E = NA_real_,
                            leverage_region = NA_character_, leverage_share = NA_real_))
  q <- c_a / D; E <- -sum(q[q > 0] * log(q[q > 0])) / log(sum(ok))
  tibble(n_regions = sum(ok), D = D, E = E,
         leverage_region = regions[ok][which.max(q)], leverage_share = max(q))
}
DE <- bind_rows(lapply(c(ST = "js_ST", SP = "js_SP", SG = "js_SG"), function(col)
  contrasts %>% group_by(pheno) %>%
    group_modify(~ de_one(.x[[col]], .x$rho, .x$target)) %>% ungroup() %>%
    mutate(pair = col)), .id = NULL) %>%
  mutate(pair = recode(pair, js_ST = "S-RAILS vs T-RAILS", js_SP = "S-RAILS vs P-RAILS",
                       js_SG = "S-RAILS vs G-RAILS"),
         leverage = ifelse(!is.na(leverage_share) & leverage_share >= 0.5, leverage_region, "spread"))
write_excel_csv(DE, out_file("DE"))

## Sign consistency of the transport contrast within each region (descriptive)
sign_cons <- contrasts %>% filter(is.finite(delta_out_within)) %>%
  group_by(target) %>%
  summarise(n_phenotypes = n(),
            pct_T_above_S = 100 * mean(delta_out_within > 0),
            median_log2_ratio_TS = median(log2_ratio_TS, na.rm = TRUE),
            pct_abs_rel_over_10 = 100 * mean(abs(rel_out_within) > 0.10, na.rm = TRUE),
            binom_p = binom.test(sum(delta_out_within > 0), n())$p.value,   # descriptive only
            .groups = "drop")
write_excel_csv(sign_cons, out_file("sign_consistency"))
message("\n--- Sign consistency of T-RAILS - S-RAILS across phenotypes ---"); print(sign_cons)

########################################################################
## PART A (cont.) — Concept shift: does Y | X differ in vs out of region a?
## Unweighted binomial fits on cells (main effects of the six covariates):
##   proxy = sum_c pop_a(c) [m_out(c) - m_in(c)] / N^a
##   LRT   = region-in/out x covariate interactions (also the effect-
##           modification test behind Part B)
########################################################################

X6 <- model.matrix(as.formula(paste("~", paste(NAMES_6, collapse = "+"))), U6)
fit_bin <- function(X, y, n) {
  ok <- n > 0
  f <- suppressWarnings(glm.fit(X[ok, , drop = FALSE], y[ok] / n[ok], weights = n[ok], family = binomial()))
  b <- f$coefficients; b[is.na(b)] <- 0
  list(coef = b, deviance = f$deviance, rank = f$rank)
}
if (RUN_CONCEPT_SHIFT) {
  message("\nConcept shift (", length(disease_cols), " phenotypes x 4 regions) ...")
  cs_list <- list()
  for (a in REGION_LEVELS) {
    ci <- src_counts(a, "within"); co <- src_counts(a, "out"); pa <- pop6[, a]
    rin <- ci$n > 0; rout <- co$n > 0
    Xp <- rbind(cbind(X6[rin, ], inr = 1), cbind(X6[rout, ], inr = 0))
    Xi <- cbind(Xp, Xp[, 2:ncol(X6)] * Xp[, "inr"])
    np <- c(ci$n[rin], co$n[rout])
    for (j in seq_along(disease_cols)) {
      yi <- ci$Y[, j]; yo <- co$Y[, j]
      if (sum(yi) < SUPPRESS_MIN || sum(yo) < SUPPRESS_MIN) next
      out <- tryCatch({
        m_in  <- fit_bin(X6, yi, ci$n); m_out <- fit_bin(X6, yo, co$n)
        p_in  <- plogis(X6 %*% m_in$coef); p_out <- plogis(X6 %*% m_out$coef)
        yp <- c(yi[rin], yo[rout])
        g0 <- fit_bin(Xp, yp, np); g1 <- fit_bin(Xi, yp, np)
        lrt <- g0$deviance - g1$deviance; df <- g1$rank - g0$rank
        tibble(target = a, pheno = disease_cols[j],
               std_within = sum(pa * p_in) / sum(pa), std_out = sum(pa * p_out) / sum(pa),
               concept_shift = sum(pa * (p_out - p_in)) / sum(pa),
               lrt = lrt, lrt_df = df, lrt_p = pchisq(lrt, df, lower.tail = FALSE))
      }, error = function(e) NULL)
      cs_list[[length(cs_list) + 1]] <- out
    }
  }
  concept_shift <- bind_rows(cs_list)
  write_excel_csv(concept_shift, out_file("concept_shift"))
}

########################################################################
## PART B — Associations with the demographic covariates
## Weighted logistic regression on cells (exact for cell-constant weights):
##   point estimate: glm.fit(X, cases/n, weights = w * n)
##   sandwich SE:    B = X' diag(w n p (1-p)) X,
##                   M = X' diag(w^2 [cases (1-p)^2 + (n - cases) p^2]) X,
##                   V = B^-1 M B^-1   (linearized; ignores calibration)
## CONDITIONAL: phenotype ~ all six covariates; MARGINAL: phenotype ~ one of
## sex / age group / race-ethnicity. Reference levels in REF_LEVELS.
########################################################################

if (RUN_BETA) {
  U6b <- U6
  for (v in names(REF_LEVELS)) U6b[[v]] <- relevel(U6b[[v]], ref = REF_LEVELS[[v]])
  DESIGNS <- list(
    conditional = model.matrix(as.formula(paste("~", paste(NAMES_6, collapse = "+"))), U6b),
    marginal_sex      = model.matrix(~ sex, U6b),
    marginal_agegroup = model.matrix(~ agegroup, U6b),
    marginal_race_eth = model.matrix(~ race_eth, U6b))
  ## Terms reported and the exposure variable / level each belongs to
  TERMS <- bind_rows(lapply(c("sex", "agegroup", "race_eth"), function(v) {
    lv <- setdiff(levels(U6b[[v]]), REF_LEVELS[[v]])
    tibble(variable = v, level = lv, term = paste0(v, lv))
  }))
  level_of <- lapply(setNames(names(REF_LEVELS), names(REF_LEVELS)), function(v) as.character(U6b[[v]]))

  fit_wlogit <- function(X, y, n, w) {
    ok <- n > 0 & w > 0
    X <- X[ok, , drop = FALSE]; y <- y[ok]; n <- n[ok]; w <- w[ok]
    f <- suppressWarnings(glm.fit(X, y / n, weights = w * n, family = binomial()))
    b <- f$coefficients; keep <- !is.na(b)
    if (!f$converged || sum(keep) == 0) return(NULL)
    Xk <- X[, keep, drop = FALSE]; p <- f$fitted.values
    B <- crossprod(Xk, Xk * (w * n * p * (1 - p)))
    M <- crossprod(Xk, Xk * (w^2 * (y * (1 - p)^2 + (n - y) * p^2)))
    Bi <- tryCatch(solve(B), error = function(e) NULL)
    if (is.null(Bi)) return(NULL)
    se <- rep(NA_real_, length(b)); se[keep] <- sqrt(pmax(diag(Bi %*% M %*% Bi), 0))
    list(beta = setNames(b, colnames(X)), se = setNames(se, colnames(X)))
  }

  combos <- bind_rows(
    expand_grid(target = REGION_LEVELS, source = SOURCES, method = c("unweighted", "rails")),
    tibble(target = REGION_LEVELS, source = "within", method = "grails"))
  message("\nPart B: ", nrow(combos), " (region x source x weighting) combinations x up to ",
          length(disease_cols), " phenotypes x ", length(DESIGNS), " models ...")
  beta_list <- list(); t0 <- tic()
  for (k in seq_len(nrow(combos))) {
    a <- combos$target[k]; s <- combos$source[k]; m <- combos$method[k]
    cnt <- src_counts(a, s)
    w <- if (m == "unweighted") as.numeric(cnt$n > 0) else cell_weights(a, s, m)
    if (all(is.na(w))) next
    cases_tot <- colSums(cnt$Y)
    ## cases per exposure level in this source (disclosure / stability rule)
    lvl_cases <- lapply(names(level_of), function(v) rowsum(cnt$Y, level_of[[v]]))
    names(lvl_cases) <- names(level_of)
    for (j in which(cases_tot >= MIN_CASES_BETA)) {
      for (mod in names(DESIGNS)) {
        fr <- fit_wlogit(DESIGNS[[mod]], cnt$Y[, j], cnt$n, w)
        if (is.null(fr)) next
        tt <- TERMS %>% filter(term %in% names(fr$beta))
        if (nrow(tt) == 0) next
        ok_lvl <- mapply(function(v, l) {
          lc <- lvl_cases[[v]][, j]
          isTRUE(lc[l] >= SUPPRESS_MIN) && isTRUE(lc[REF_LEVELS[[v]]] >= SUPPRESS_MIN)
        }, tt$variable, tt$level)
        b  <- fr$beta[tt$term]; se <- fr$se[tt$term]
        bad <- !ok_lvl | !is.finite(b) | abs(b) > BETA_MAX_ABS
        beta_list[[length(beta_list) + 1]] <- tibble(
          target = a, source = s, method = m, model = ifelse(mod == "conditional", "conditional", "marginal"),
          pheno = disease_cols[j], variable = tt$variable, term = tt$term,
          beta = ifelse(bad, NA_real_, b), se = ifelse(bad, NA_real_, se))
      }
    }
    message(format(Sys.time(), "%H:%M:%S"), "  [", k, "/", nrow(combos), "] ", a, " / ", s, " / ", m,
            "  (", round(tic() - t0), " s)")
  }
  beta <- bind_rows(beta_list)
  write_excel_csv(beta, out_file("beta"))

  ## Transport contrasts in beta (RAILS weights), with an approximate z.
  ## within and out are disjoint samples, so var(T - S) ~ se_T^2 + se_S^2;
  ## all overlaps within, so its z is only indicative.
  beta_contrasts <- beta %>%
    filter(method %in% c("rails", "unweighted")) %>%
    select(target, method, model, pheno, variable, term, source, beta, se) %>%
    pivot_wider(names_from = source, values_from = c(beta, se))
  for (col in c(outer(c("beta_", "se_"), SOURCES, paste0)))      # a source that failed everywhere
    if (!col %in% names(beta_contrasts)) beta_contrasts[[col]] <- NA_real_
  beta_contrasts <- beta_contrasts %>%
    mutate(dbeta_out_within = beta_out - beta_within,
           z_out_within     = dbeta_out_within / sqrt(se_out^2 + se_within^2),
           dbeta_all_within = beta_all - beta_within,
           z_all_within     = dbeta_all_within / sqrt(se_all^2 + se_within^2)) %>%
    left_join(contrasts %>% select(target, pheno, logit_delta_TS), by = c("target", "pheno"))
  write_excel_csv(beta_contrasts, out_file("beta_contrasts"))

  ## Summary: is the association more transportable than the prevalence?
  ## Both on the log-odds scale: |logit T - logit S| vs |beta_T - beta_S|.
  beta_summary <- bind_rows(
    beta_contrasts %>% filter(method == "rails", is.finite(dbeta_out_within)) %>%
      group_by(model, variable, term, target) %>%
      summarise(n_phenotypes = n(), median_abs_delta = median(abs(dbeta_out_within)),
                pct_abs_z_over_1.96 = 100 * mean(abs(z_out_within) > 1.96, na.rm = TRUE),
                .groups = "drop"),
    contrasts %>% filter(is.finite(logit_delta_TS)) %>% group_by(target) %>%
      summarise(n_phenotypes = n(), median_abs_delta = median(abs(logit_delta_TS)), .groups = "drop") %>%
      mutate(model = "prevalence (logit)", variable = NA_character_, term = NA_character_,
             pct_abs_z_over_1.96 = NA_real_))
  write_excel_csv(beta_summary, out_file("beta_summary"))
  message("\n--- Median |T - S| on the log-odds scale: prevalence vs associations (RAILS) ---")
  print(beta_summary %>% group_by(model) %>%
          summarise(median_abs_delta = median(median_abs_delta), .groups = "drop"))
}

########################################################################
## Figures (../graph/transport2_fig*.png) — same visual system as T01
########################################################################

INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; SURFACE <- "#fcfcfb"
SRC_COLS    <- c(within = "#2a78d6", all = "#eb6834", out = "#1baf7a")
REGION_COLS <- c(Northeast = "#2a78d6", Midwest = "#eb6834", South = "#1baf7a", West = "#eda100")
REGION_SHAPES <- c(Northeast = 16, Midwest = 17, South = 15, West = 18)   # secondary encoding
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
save_fig <- function(g, name, w, h) {
  f <- file.path(GRAPH_DIR, paste0("transport2_", name, ".png"))
  ggsave(f, g, width = w, height = h, dpi = 200, bg = SURFACE); message("Wrote ", f)
}
CAP <- "S-, P-, T-RAILS = RAILS on the within-region (U_A^a), all-region (U_A) and out-of-region (U_A^{-a}) source; G-RAILS = manuscript's global weights. No external benchmark."
reg_f <- function(x) factor(x, REGION_LEVELS)

## fig1: distribution of log2(T/S) and log2(P/S) across phenotypes, by region
f1 <- contrasts %>%
  select(target, pheno, `T-RAILS / S-RAILS` = log2_ratio_TS, `P-RAILS / S-RAILS` = log2_ratio_PS) %>%
  pivot_longer(-c(target, pheno), names_to = "contrast", values_to = "l2") %>%
  filter(is.finite(l2)) %>% mutate(target = reg_f(target),
                                   contrast = factor(contrast, c("T-RAILS / S-RAILS", "P-RAILS / S-RAILS")))
g1 <- ggplot(f1, aes(x = target, y = l2, colour = contrast)) +
  geom_hline(yintercept = 0, colour = INK2, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 0.7, alpha = 0.35) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.5, fill = NA, outlier.shape = NA, linewidth = 0.5) +
  scale_colour_manual(values = c(`T-RAILS / S-RAILS` = SRC_COLS[["out"]], `P-RAILS / S-RAILS` = SRC_COLS[["all"]]),
                      name = NULL) +
  scale_y_continuous(breaks = -3:3, labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8")) +
  labs(x = "Target region", y = "Ratio to S-RAILS (log2 scale)",
       title = "How much does the source sample change each phenotype's prevalence?",
       subtitle = "One dot per phenotype; box = median and interquartile range across phenotypes.",
       caption = CAP) +
  theme_transport()
save_fig(g1, "fig1_prevalence_ratio_distribution", 10, 5.5)

## fig2: sign consistency of the transport contrast
g2 <- ggplot(sign_cons %>% mutate(target = reg_f(target)), aes(x = target, y = pct_T_above_S)) +
  geom_hline(yintercept = 50, colour = INK2, linewidth = 0.4) +
  geom_col(width = 0.55, fill = SRC_COLS[["out"]]) +
  geom_text(aes(label = sprintf("%.0f%%\n(n = %d)", pct_T_above_S, n_phenotypes)), vjust = -0.3,
            size = 3, colour = INK) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  labs(x = "Target region", y = "% of phenotypes with T-RAILS > S-RAILS",
       title = "Is the direction of the transport contrast shared across phenotypes?",
       subtitle = "50% = no common direction. A share far from 50% points to a region-level factor (e.g. ascertainment), not phenotype-specific epidemiology.",
       caption = CAP) +
  theme_transport()
save_fig(g2, "fig2_sign_consistency", 8, 5)

## fig3: transport contrast vs concept-shift proxy, all phenotypes
if (RUN_CONCEPT_SHIFT && nrow(concept_shift) > 0) {
  f3 <- contrasts %>% select(target, pheno, delta_out_within) %>%
    inner_join(concept_shift %>% select(target, pheno, concept_shift, lrt_p), by = c("target", "pheno")) %>%
    filter(is.finite(delta_out_within), is.finite(concept_shift)) %>%
    mutate(target = reg_f(target), sig = ifelse(lrt_p < 0.05, "p < 0.05", "p >= 0.05"))
  lim3 <- f3 %>% group_by(target) %>%
    summarise(l = max(abs(100 * c(delta_out_within, concept_shift)), 0.05) * 1.1, .groups = "drop") %>%
    reframe(x = c(-l, l), y = c(-l, l), .by = target)
  g3 <- ggplot(f3, aes(x = 100 * concept_shift, y = 100 * delta_out_within)) +
    geom_blank(data = lim3, aes(x = x, y = y)) +
    geom_hline(yintercept = 0, colour = GRID) + geom_vline(xintercept = 0, colour = GRID) +
    geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4) +
    geom_point(aes(fill = sig), shape = 21, colour = INK, size = 1.6, stroke = 0.3, alpha = 0.8) +
    scale_fill_manual(values = c(`p < 0.05` = INK, `p >= 0.05` = SURFACE), limits = c("p < 0.05", "p >= 0.05"),
                      name = "In/out covariate effects differ (LRT)") +
    facet_wrap(~ target, scales = "free", nrow = 1) +
    labs(x = "Concept-shift proxy: mean over target of [m_out(X) - m_within(X)] (pp)",
         y = "T-RAILS - S-RAILS (pp)",
         title = "Is the transport contrast explained by concept shift, phenotype by phenotype?",
         subtitle = "One dot per phenotype. On the diagonal: the contrast is what the regional difference in outcome given X predicts.",
         caption = paste(CAP, "Both axes come from the same AoU data, so agreement is partly expected by construction.")) +
    theme_transport() + theme(aspect.ratio = 1)
  save_fig(g3, "fig3_contrast_vs_concept_shift", 12, 4.8)
}

## fig4: D vs E (manuscript metrics, unit = region) for S-T, S-P and S-G
f4 <- DE %>% filter(is.finite(D), is.finite(E), D > 0) %>%
  mutate(pair = factor(pair, c("S-RAILS vs T-RAILS", "S-RAILS vs P-RAILS", "S-RAILS vs G-RAILS")),
         leverage = factor(leverage, c(REGION_LEVELS, "spread")))
g4 <- ggplot(f4, aes(x = D, y = E, colour = leverage, shape = leverage)) +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_x_log10() +
  scale_colour_manual(values = c(REGION_COLS, spread = "#8a8984"), name = "Leverage region (share >= 50%)") +
  scale_shape_manual(values = c(REGION_SHAPES, spread = 1), name = "Leverage region (share >= 50%)") +
  facet_wrap(~ pair, nrow = 1) +
  labs(x = expression(D[j] ~ "(population-weighted Jensen-Shannon divergence, log scale)"),
       y = expression(E[j] ~ "(dispersion across the 4 regions)"),
       title = "How large is each source contrast, and is it concentrated in one region?",
       subtitle = "One dot per phenotype. E: 1 = spread evenly over the 4 regions, 0 = one region carries all of D (coarse with L = 4).",
       caption = CAP) +
  theme_transport()
save_fig(g4, "fig4_D_vs_E", 12, 5)

if (RUN_BETA && nrow(beta) > 0) {
  term_lab <- function(t) sub("^(sex|agegroup|race_eth)", "", t)
  var_lab  <- c(sex = "Sex (ref. Female)", agegroup = "Age (ref. 45-64)", race_eth = "Race/ethnicity (ref. NH White)")

  ## fig5: beta(T) vs beta(S), conditional model, by term
  f5 <- beta_contrasts %>% filter(method == "rails", model == "conditional",
                                  is.finite(beta_within), is.finite(beta_out)) %>%
    mutate(target = reg_f(target), panel = paste0(var_lab[variable], ": ", term_lab(term)))
  g5 <- ggplot(f5, aes(x = beta_within, y = beta_out, colour = target, shape = target)) +
    geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4) +
    geom_point(size = 1.2, alpha = 0.6) +
    scale_colour_manual(values = REGION_COLS, name = "Target region") +
    scale_shape_manual(values = REGION_SHAPES, name = "Target region") +
    facet_wrap(~ panel, scales = "free", ncol = 3) +
    labs(x = "beta under S-RAILS (within-region source)", y = "beta under T-RAILS (out-of-region source)",
         title = "Do adjusted demographic associations transport?",
         subtitle = "Conditional log-odds ratios (phenotype ~ all six covariates). One dot per phenotype x region; on the diagonal = same association.",
         caption = CAP) +
    theme_transport()
  save_fig(g5, "fig5_beta_T_vs_S", 11, 9)

  ## fig6: transport discrepancy on the log-odds scale: prevalence vs marginal vs conditional beta
  f6 <- bind_rows(
    contrasts %>% filter(is.finite(logit_delta_TS)) %>%
      transmute(target, quantity = "Prevalence: |logit T - logit S|", abs_delta = abs(logit_delta_TS)),
    beta_contrasts %>% filter(method == "rails", is.finite(dbeta_out_within)) %>%
      transmute(target, quantity = ifelse(model == "marginal", "Marginal beta: |beta_T - beta_S|",
                                          "Conditional beta: |beta_T - beta_S|"),
                abs_delta = abs(dbeta_out_within))) %>%
    mutate(target = reg_f(target),
           quantity = factor(quantity, c("Prevalence: |logit T - logit S|", "Marginal beta: |beta_T - beta_S|",
                                         "Conditional beta: |beta_T - beta_S|")))
  g6 <- ggplot(f6, aes(x = target, y = abs_delta, colour = quantity)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.size = 0.4, outlier.alpha = 0.3) +
    scale_colour_manual(values = c("#52514e", "#eb6834", "#2a78d6"), name = NULL) +
    scale_y_sqrt() +
    labs(x = "Target region", y = "Absolute T - S difference, log-odds scale (square-root axis)",
         title = "Do associations transport better than prevalences?",
         subtitle = "All three on the log-odds scale. Prevalence: one value per phenotype; beta: one per phenotype x demographic term.",
         caption = CAP) +
    theme_transport() + guides(colour = guide_legend(nrow = 1))
  save_fig(g6, "fig6_prevalence_vs_beta_transport", 10, 5.5)

  ## fig7: share of significant beta contrasts, by term x region, marginal and conditional
  f7 <- beta_summary %>% filter(model != "prevalence (logit)") %>%
    mutate(target = reg_f(target), term = paste0(var_lab[variable], ": ", term_lab(term)),
           txt_col = ifelse(coalesce(pct_abs_z_over_1.96, 0) > 60, "#ffffff", INK))
  g7 <- ggplot(f7, aes(x = target, y = term, fill = pct_abs_z_over_1.96)) +
    geom_tile(colour = SURFACE, linewidth = 1) +
    geom_text(aes(label = ifelse(is.finite(pct_abs_z_over_1.96), sprintf("%.0f%%", pct_abs_z_over_1.96), "–"),
                  colour = txt_col), size = 2.8) +
    scale_colour_identity() +
    scale_fill_gradient(low = "#f0efec", high = "#184f95", limits = c(0, 100), name = "% |z| > 1.96") +
    facet_wrap(~ model) +
    labs(x = "Target region", y = NULL,
         title = "How often does T-RAILS give a different association than S-RAILS?",
         subtitle = "Share of phenotypes with |beta_T - beta_S| / SE > 1.96. With naive SEs about 5% is expected by chance if associations transport.",
         caption = CAP) +
    theme_transport() + theme(panel.grid.major = element_blank(), legend.position = "right")
  save_fig(g7, "fig7_beta_contrast_significance", 11, 5.5)
}

message("\nDone. Tables: ", DATA_DIR, "/transport2_*.csv | figures: ", GRAPH_DIR, "/transport2_fig*.png",
        "\n(transport2_cells.csv contains raw cell counts — keep it inside the Workbench.)")
