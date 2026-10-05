#### VEHSS_comparison_cataract.R — Cataract, adults 40+: AoU (unweighted /
#### Global RAILS / Subgroup RAILS by region) versus the CDC Vision & Eye Health
#### Surveillance System (VEHSS) state estimates.
####
#### Same pipeline as VEHSS_comparison.R (AMD), phenotype swapped: every stage
#### re-runs on each execution (SKIP_IF_CACHED <- TRUE only to reuse stage
#### outputs); the VEHSS reference tables are NOT downloaded here — run
#### VEHSS_data_eye.R where internet access exists; it writes VEHSS_CATARACT_overall.csv
#### and VEHSS_CATARACT_combinations.csv, which this script reads from its own
#### directory (or ../data/). Every output is prefixed vehss_cataract_ so nothing
#### from the AMD run is overwritten (the raw ACS download cache is shared).
#### Stages: 1 PUMS 40+ by state (Census API) | 2 AoU 40+ cohort + the ICD-
#### defined phenotype (Workbench 2.0 / BigQuery) | 3 Global RAILS + Subgroup
#### RAILS by region, Global RAILS fallback in a region that fails | 4 state /
#### region / national comparison, D_j / E_j per Age x Sex x Race cell, maps,
#### scatter, dot plot, D-vs-E. Covariates, age cuts (40-64 / 65-84 / 85+),
#### calibration and disclosure rule (SUPPRESS_MIN) are identical to AMD.
####
#### Reference: VEHSS claims source "Medicare Fee for Service and Medicare
#### Advantage", indicator "Annual prevalence of diagnosed cataracts" /
#### "All Cataracts" (Diagnosed Cataract, Any Type - Any Stage), risk factor
#### "All patients" (the diabetes / hypertension splits are NOT used), year
#### 2022 (the portal id YR19 in the given URL is 2022, not 2019), crude
#### prevalence, by state + US. Differences from the AMD reference:
####   * DIAGNOSED prevalence from claims (one code suffices), not a modeled
####     population estimate — closer in kind to the AoU EHR definition;
####   * DENOMINATOR = Medicare beneficiaries, not residents (implied 40+
####     population 52.8M vs 162M in ACS): 65+ is near-universal coverage,
####     but 40-64 is disability / ESRD beneficiaries only, so the 40-64
####     cells are not comparable to AoU's general 40-64 population;
####   * ages offered: all | 0-17 | 18-39 | 18+ | 40-64 | 40+ | 65-84 | 65+ |
####     85+ (the four AMD codes are used); sexes as AMD;
####   * races: White NH | Black NH | Hispanic | Asian | North American Native
####     | Other, so the portal's "Other" EXCLUDES Asian and AIAN. The
####     comparable "Other" in VEHSS_CATARACT_*.csv is DERIVED (Other + Asian
####     + AIAN, or All - White - Black - Hispanic where a component is
####     suppressed; no CI); the three portal groups are kept as extra levels.
####
#### AoU case definition (ICD-9-CM / ICD-10-CM SOURCE codes on >= MIN_CODE_DATES
#### distinct dates, prefix match; searched in condition_occurrence AND
#### observation, because OMOP routes status codes such as Z96.1 "presence of
#### intraocular lens" to the observation table): VEHSS crosswalk category 6
#### "Cataracts", ALL subgroups, as the reference's "All Cataracts" is —
####   6.1 senile 366.1x / H25; 6.2 non-congenital 366.2-366.4, 366.8, 366.9 /
####   H26.1-H26.3, H26.8, H26.9, H28; 6.3 congenital 366.0x, 743.3x / H26.0x,
####   Q12.0/.2/.3/.4/.8; 6.4 posterior capsular opacity 366.5x / H26.4x;
####   6.5 pseudophakia V43.1 / Z96.1; 6.6 aphakia & other lens disorders
####   379.3x / H27. Subgroups 6.3 and 6.5-6.6 are switchable below.
####
#### Requires Sub_AoU_Fun.R AND AoU_Fun.R next to this script (stage 3).
#### Inputs : VEHSS_CATARACT_overall.csv, VEHSS_CATARACT_combinations.csv (from
####          VEHSS_data_eye.R), Census API key (optional, env CENSUS_API_KEY)
####          for stage 1, AoU CDR for stage 2.
#### Outputs: ../data/vehss_cataract_*.csv, ../graph/vehss_cataract_*.png
####
#### AoU DISCLOSURE POLICY: any state with fewer than SUPPRESS_MIN cases is
#### suppressed (hatched on maps, NA in tables). Maps are contiguous-US only.

suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(Matrix)
})
## `maps` supplies the state polygons; `ggrepel` / `ggpattern` are optional.
if (!requireNamespace("maps", quietly = TRUE)) install.packages("maps")

########################################################################
## Settings
########################################################################

SKIP_IF_CACHED     <- FALSE   # FALSE (default) = re-run every stage; TRUE = reuse existing stage outputs
REDOWNLOAD_PUMS_RAW <- FALSE  # TRUE = fetch the raw ACS PUMS records from the Census API again
                              # (the raw download is data, not a processing step; the recoding,
                              #  the 40+ restriction and everything after it always re-run)
AGE_MIN        <- 40      # VEHSS "40 years and older"
## Age-group cuts INSIDE the 40+ window (left-closed). Default = VEHSS strata.
## A finer alternative that keeps the pipeline's upper cuts: c(40, 50, 65, 75, Inf)
AGE_BREAKS     <- c(40, 65, 85, Inf)
MIN_CODE_DATES <- 2       # distinct code dates required to be a case
SUPPRESS_MIN   <- 20      # AoU small-cell threshold (cases per state)
ALPHA          <- 0.05    # forward LRT threshold for the three-way selection

## VEHSS "All Cataracts" (R6_ALL) is the whole crosswalk category 6, INCLUDING
## pseudophakia (Sub_6.5), aphakia / other lens disorders (Sub_6.6) and
## congenital cataract (Sub_6.3). Keep both TRUE to match the reference.
CATARACT_INCLUDE_PSEUDOPHAKIA_APHAKIA <- TRUE
CATARACT_INCLUDE_CONGENITAL           <- TRUE
PHENO_ID    <- "cataract"
PHENO_LABEL <- "Cataract"
PHENO_ICD9  <- c("366.1", "366.2", "366.3", "366.4", "366.5", "366.8", "366.9",
                 if (CATARACT_INCLUDE_CONGENITAL) c("366.0", "743.3"),
                 if (CATARACT_INCLUDE_PSEUDOPHAKIA_APHAKIA) c("V43.1", "379.3"))
PHENO_ICD10 <- c("H25", "H26.1", "H26.2", "H26.3", "H26.4", "H26.8", "H26.9", "H28",
                 if (CATARACT_INCLUDE_CONGENITAL) c("H26.0", "Q12.0", "Q12.2", "Q12.3", "Q12.4", "Q12.8"),
                 if (CATARACT_INCLUDE_PSEUDOPHAKIA_APHAKIA) c("Z96.1", "H27"))

## ---- Case-definition switch ---------------------------------------------
## "phecodeX": the PhecodeX phecode (exact ICD-9-CM / ICD-10-CM codes, child
##             phecodes included) — the same phecode-based phenotyping as the
##             manuscript's phecode categories (06_Phecode_Categories.R)
## "vehss":    the VEHSS-crosswalk-style prefix lists above
## Differences for cataract: PhecodeX SO_371 INCLUDES diabetic cataract
## E08-E13.36, retained cataract fragments after surgery 998.82 / H59.02x and
## cataract-extraction status V45.61 / Z98.4x; it EXCLUDES pseudophakia Z96.1
## (the ICD-9 V43.1 is included), aphakia 379.31 / H27 (PhecodeX SO_372.1),
## congenital cataract 743.31-743.36 / Q12.2-Q12.8 (CM_751.213), H26.11x and
## 366.21, which the VEHSS "All Cataracts" list carries.
PHENO_DEF <- "phecodeX"
## PhecodeX SO_371 (with its child phecodes), from phecodeX_R_map.csv + phecodeX_R_rollup_map.csv
## (github.com/PheWAS/PhecodeX, downloaded 2026-09-18): exact ICD codes, 48 ICD-9-CM + 141 ICD-10-CM
PHECODEX_CATARACT <- list(phecode = "SO_371",
  icd9 = c(
    "366", "366.0", "366.00", "366.01", "366.02", "366.03", "366.04", "366.09", "366.1", "366.10",
    "366.11", "366.12", "366.13", "366.14", "366.15", "366.16", "366.17", "366.18", "366.19",
    "366.2", "366.20", "366.22", "366.23", "366.3", "366.30", "366.31", "366.32", "366.33", "366.34",
    "366.4", "366.41", "366.42", "366.43", "366.44", "366.45", "366.46", "366.5", "366.50", "366.51",
    "366.52", "366.53", "366.8", "366.9", "743.3", "743.30", "998.82", "V43.1", "V45.61"),
  icd10 = c(
    "E08.36", "E09.36", "E10.36", "E11.36", "E13.36", "H25", "H25.0", "H25.01", "H25.011", "H25.012",
    "H25.013", "H25.019", "H25.03", "H25.031", "H25.032", "H25.033", "H25.039", "H25.04", "H25.041",
    "H25.042", "H25.043", "H25.049", "H25.09", "H25.091", "H25.092", "H25.093", "H25.099", "H25.1",
    "H25.10", "H25.11", "H25.12", "H25.13", "H25.2", "H25.20", "H25.21", "H25.22", "H25.23", "H25.8",
    "H25.81", "H25.811", "H25.812", "H25.813", "H25.819", "H25.89", "H25.9", "H26", "H26.0",
    "H26.00", "H26.001", "H26.002", "H26.003", "H26.009", "H26.01", "H26.011", "H26.012", "H26.013",
    "H26.019", "H26.03", "H26.031", "H26.032", "H26.033", "H26.039", "H26.04", "H26.041", "H26.042",
    "H26.043", "H26.049", "H26.05", "H26.051", "H26.052", "H26.053", "H26.059", "H26.06", "H26.061",
    "H26.062", "H26.063", "H26.069", "H26.09", "H26.1", "H26.10", "H26.101", "H26.102", "H26.103",
    "H26.109", "H26.12", "H26.121", "H26.122", "H26.123", "H26.129", "H26.13", "H26.131", "H26.132",
    "H26.133", "H26.139", "H26.2", "H26.20", "H26.21", "H26.211", "H26.212", "H26.213", "H26.219",
    "H26.22", "H26.221", "H26.222", "H26.223", "H26.229", "H26.23", "H26.231", "H26.232", "H26.233",
    "H26.239", "H26.3", "H26.30", "H26.31", "H26.32", "H26.33", "H26.4", "H26.40", "H26.41",
    "H26.411", "H26.412", "H26.413", "H26.419", "H26.49", "H26.491", "H26.492", "H26.493", "H26.499",
    "H26.8", "H26.9", "H28", "H59.02", "H59.021", "H59.022", "H59.023", "H59.029", "Q12.0", "Z98.4",
    "Z98.41", "Z98.42", "Z98.49"))

PHECODEX_CODES <- c(PHECODEX_CATARACT$icd9, PHECODEX_CATARACT$icd10)
## SQL condition on an ICD-code column for the active case definition:
##   phecodeX -> exact membership in the phecode's code list
##   vehss    -> the prefix regex built from the lists above
code_cond <- function(col, rx) {
  if (PHENO_DEF == "phecodeX")
    paste0(col, " IN (", paste0("'", PHECODEX_CODES, "'", collapse = ", "), ")")
  else paste0("REGEXP_CONTAINS(", col, ", r'", rx, "')")
}

DATA_DIR  <- "../data"
GRAPH_DIR <- "../graph"
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(GRAPH_DIR, showWarnings = FALSE, recursive = TRUE)

## VEHSS reference tables, produced by VEHSS_data.R (looked up in this
## directory first, then in ../data/)
VEHSS_OVERALL_FILE <- "VEHSS_CATARACT_overall.csv"        # 40+ / both sexes / all races
VEHSS_COMBO_FILE   <- "VEHSS_CATARACT_combinations.csv"   # every other Age x Sex x Race selection

## Covariate-combination analysis (Age x Sex x Race, VEHSS categories)
N_SELECT   <- 4        # combinations per selection criterion (most populated, high D, low E, ...)
MIN_STATES <- 10       # a combination needs >= this many un-suppressed states for D / E
MAP_SETS   <- c("pop", "highD", "lowE")   # which selection sets get their own map figures
                                          # (any of "pop", "highD", "lowD", "highE", "lowE")
## Combinations compared: the overall (40+ / both sexes / all races) plus every
## FULLY specified cell (a specific age band x sex x race group). Partial
## combinations that leave one or two dimensions at "all" (e.g. 40-64 / both
## sexes / White, or 40+ / Female / all races) are dropped unless this is TRUE.
ALLOW_PARTIAL_COMBOS <- FALSE
## Target 40+ population the calibrated weights are scaled to. NA = the ACS
## 2022 40+ total of the states in state_region_map (sum of PWGTP over ALL
## 40+ PUMS records of those states, before dropping rows with a missing
## covariate). Set a number to force another anchor, e.g. the VEHSS-implied
## 2019 total (vehss_us$pop_vehss) or a Census population estimate.
TARGET_POP <- NA_real_
CENSUS_KEY <- Sys.getenv("CENSUS_API_KEY", unset = "")
AGE_REFERENCE_DATE <- as.Date(Sys.getenv("AGE_REFERENCE_DATE", unset = "2024-08-01"))

## Stage outputs (all `vehss_`-prefixed so nothing in the main pipeline is touched)
F_PUMS_RAW <- file.path(DATA_DIR, "vehss_PUMS_2022_state_raw.csv")   # raw API download, all ages
F_PUMS     <- file.path(DATA_DIR, "vehss_cataract_pums_2022_bystate_40plus.csv")
F_AOU      <- file.path(DATA_DIR, "vehss_cataract_aou_40plus.csv")
F_WEIGHTS  <- file.path(DATA_DIR, "vehss_cataract_rails_weights_40plus.csv")
F_CELLS_G  <- file.path(DATA_DIR, "vehss_cataract_rails_cells_global_40plus.csv")
F_CELLS_S  <- file.path(DATA_DIR, "vehss_cataract_rails_cells_region_40plus.csv")
F_STATE    <- file.path(DATA_DIR, "vehss_cataract_state_comparison.csv")
F_REGION   <- file.path(DATA_DIR, "vehss_cataract_region_comparison.csv")
F_NATIONAL <- file.path(DATA_DIR, "vehss_cataract_national_comparison.csv")
F_SUMMARY  <- file.path(DATA_DIR, "vehss_cataract_method_summary.csv")

## RAILS function files — located NOW, before the long stages, so a missing
## file cannot kill the run after the BigQuery stage has already been paid
## for. Sub_AoU_Fun.R and AoU_Fun.R must sit in the same directory (the
## working directory, its parent, ../data, the repo folders, or anywhere below
## the home / workspace directory).
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
    stop(name, " not found. Upload it (together with AoU_Fun.R, from Application/Global RAILS/",
         "RAILS Procedure/) into ", normalizePath(getwd()), " and re-run.")
  normalizePath(f[1])
}
SUB_FUN_FILE <- find_fun_file("Sub_AoU_Fun.R")
if (!file.exists(file.path(dirname(SUB_FUN_FILE), "AoU_Fun.R")))
  stop("AoU_Fun.R must sit next to ", SUB_FUN_FILE, " (Sub_AoU_Fun.R sources it) — upload it there.")
message("RAILS functions: ", SUB_FUN_FILE, " (+ AoU_Fun.R)")

########################################################################
## Shared definitions: covariates, age groups, factor levels, state -> region
########################################################################

NAMES_7 <- c("agegroup", "sex", "edu", "homeown", "income", "race_eth", "region")
NAMES_6 <- setdiff(NAMES_7, "region")   # within-region model for Subgroup RAILS

## Labels "40-64", "65-84", "85+" from AGE_BREAKS
AGE_LABELS <- {
  b <- AGE_BREAKS
  k <- length(b) - 1
  c(paste0(b[seq_len(k - 1)], "-", b[seq_len(k - 1) + 1] - 1), paste0(b[k], "+"))
}
make_agegroup <- function(age) {
  as.character(cut(age, breaks = AGE_BREAKS, right = FALSE, labels = AGE_LABELS))
}

## State -> region map: the OFFICIAL Census regions, identical to the PUMS `REGION` variable
## (01_PUMS_Prep.R) and to the maps of 02_AoU_Prep.R / 04_Global_RAILS.R / 09_Sub_RAILS_*.R:
## Delaware, Maryland, DC and Oklahoma are in the South. Applied to BOTH the AoU and the
## by-state PUMS side here (and to the VEHSS state benchmark), so region is defined
## identically in all three.
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

## Factor levels shared by AoU and PUMS (agegroup uses the 40+ cuts above)
harmonize_factors <- function(df) {
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
    ) %>%
    na.omit()
}

## State FIPS -> USPS (50 states + DC); drives the PUMS per-state download
fips2usps <- c(
  "01"="AL","02"="AK","04"="AZ","05"="AR","06"="CA","08"="CO","09"="CT",
  "10"="DE","11"="DC","12"="FL","13"="GA","15"="HI","16"="ID","17"="IL",
  "18"="IN","19"="IA","20"="KS","21"="KY","22"="LA","23"="ME","24"="MD",
  "25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT","31"="NE",
  "32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND",
  "39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD",
  "47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV",
  "55"="WI","56"="WY")

## A stage is skipped ONLY when reuse was explicitly requested and its output exists
stage_cached <- function(f) {
  if (SKIP_IF_CACHED && file.exists(f)) {
    message("SKIP_IF_CACHED = TRUE: reusing ", f); TRUE
  } else FALSE
}

########################################################################
## VEHSS reference tables, built SEPARATELY by VEHSS_data.R (run it on a
## machine with internet access, then upload the two csvs next to this script):
##   VEHSS_CATARACT_overall.csv       40+ / both sexes / all races, by state + US
##   VEHSS_CATARACT_combinations.csv  every other Age x Sex x Race selection
## Columns used: age, sex, race, state, prev / lb / ub (proportions), cases
## (estimated number with Cataract) and pop (implied population = cases / prev).
## VEHSS categories: age 40+ | 40-64 | 65-84 | 85+; sex Both | Female | Male;
## race All | Black | Hispanic | White | Other.
########################################################################

AGE_V  <- c("40+", "40-64", "65-84", "85+")
SEX_V  <- c("Both", "Female", "Male")
RACE_V <- c("All", "Black", "Hispanic", "White", "Other")

find_input <- function(name) {
  cands <- unique(c(name, file.path(DATA_DIR, name)))
  f <- cands[file.exists(cands)]
  if (length(f) == 0)
    stop(name, " not found (looked in ", paste(cands, collapse = ", "),
         ") — run VEHSS_data.R first and copy its outputs next to this script.")
  f[1]
}
read_vehss <- function(f) {
  readr::read_csv(f, show_col_types = FALSE,
                  col_types = readr::cols(age = "c", sex = "c", race = "c", state = "c",
                                          .default = readr::col_guess())) %>%
    transmute(age, sex, race, state,
              prev_vehss = prev, lb_vehss = lb, ub_vehss = ub,
              cases_vehss = as.numeric(cases), pop_vehss = as.numeric(pop))
}
vehss_overall <- read_vehss(find_input(VEHSS_OVERALL_FILE))
vehss_combos  <- read_vehss(find_input(VEHSS_COMBO_FILE))
vehss_all     <- bind_rows(vehss_overall, vehss_combos)
if (nrow(vehss_overall) == 0) stop("VEHSS overall table is empty.")
message("VEHSS: overall table ", nrow(vehss_overall), " rows; combinations table ",
        nrow(vehss_combos), " rows covering ",
        n_distinct(paste(vehss_combos$age, vehss_combos$sex, vehss_combos$race)),
        " Age x Sex x Race selections")

vehss_us    <- vehss_overall %>% filter(state == "US")
vehss_state <- vehss_overall %>% filter(state != "US") %>%
  mutate(region = unname(state_region_map[state]))
vehss_region <- vehss_state %>%
  filter(!is.na(region), !is.na(prev_vehss)) %>%
  group_by(region) %>%
  summarise(prev_vehss  = sum(cases_vehss) / sum(pop_vehss),
            cases_vehss = sum(cases_vehss),
            pop_vehss   = sum(pop_vehss), .groups = "drop") %>%
  mutate(lb_vehss = NA_real_, ub_vehss = NA_real_)   # modeled state CIs do not aggregate

message("VEHSS national Cataract prevalence (40+, reference VEHSS Medicare FFS+MA diagnosed (2022)): ",
        sprintf("%.2f%% (%.2f-%.2f)", 100 * vehss_us$prev_vehss,
                100 * vehss_us$lb_vehss, 100 * vehss_us$ub_vehss),
        " | implied 40+ population: ", format(round(vehss_us$pop_vehss), big.mark = ","))

########################################################################
## STAGE 1 — PUMS 2022 by state, age >= AGE_MIN
## Mirrors PUMS_2022_bystate.R (per-state ucgid calls, ST -> USPS, official
## codings) but (a) keeps only age >= AGE_MIN, (b) re-cuts agegroup with
## AGE_BREAKS, (c) derives region from state_region_map so the PUMS region is
## defined exactly like the AoU region, and (d) writes to a vehss_-prefixed file.
## PWGTP is kept raw over ALL 40+ records (covariates may be NA): that total
## is the population anchor the calibrated weights are rescaled to.
########################################################################

if (!stage_cached(F_PUMS)) {

  library(httr)
  library(jsonlite)

  ## Reuse the raw all-ages state download of PUMS_2022_bystate.R if present
  ## (same API extract; it is only read, never written)
  pipeline_raw <- file.path(DATA_DIR, "PUMS_2022_state.csv")
  raw_src <- if (file.exists(F_PUMS_RAW)) F_PUMS_RAW else if (file.exists(pipeline_raw)) pipeline_raw else NA

  GET_VARS <- paste0("PWGTP,HINCP,AGEP,RACSOR,RACAIAN,RACASN,RACBLK,RACWHT,",
                     "TEN,SEX,RACNH,HISP,SCHL,RACPI,REGION,ST")

  fetch_pums_state <- function(ucgid) {
    key_q <- if (nzchar(CENSUS_KEY)) paste0("&key=", CENSUS_KEY) else ""
    url   <- paste0("https://api.census.gov/data/2022/acs/acs1/pums?get=",
                    GET_VARS, "&ucgid=", ucgid, key_q)
    resp  <- GET(url)
    txt   <- content(resp, as = "text", encoding = "UTF-8")
    if (http_error(resp) || !startsWith(trimws(txt), "[")) {
      stop("Census API error for ucgid=", ucgid, " (HTTP ", status_code(resp), ").\n",
           substr(trimws(txt), 1, 300), call. = FALSE)
    }
    m  <- fromJSON(txt)
    df <- as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE)
    colnames(df) <- m[1, ]
    df[, !duplicated(toupper(colnames(df))), drop = FALSE]   # ST comes back twice
  }

  if (REDOWNLOAD_PUMS_RAW || is.na(raw_src)) {
    ucgids <- paste0("0400000US", names(fips2usps))
    message("Downloading PUMS 2022 by state (", length(ucgids), " calls)...")
    parts <- lapply(seq_along(ucgids), function(i) {
      message("  [", i, "/", length(ucgids), "] ", fips2usps[i])
      fetch_pums_state(ucgids[i])
    })
    dt_list <- dplyr::bind_rows(parts)
    write.csv(dt_list, F_PUMS_RAW, row.names = FALSE)
  } else {
    message("Reading cached raw PUMS: ", raw_src)
    dt_list <- read.csv(raw_src)
  }
  colnames(dt_list) <- sub("\\.\\.\\..*$", "", colnames(dt_list))
  dt_list <- dt_list[, !duplicated(toupper(colnames(dt_list)))]

  pums40 <- dt_list %>%
    rename(income = HINCP, age = AGEP, race_Asian = RACASN, race_Black = RACBLK,
           race_White = RACWHT, sex = SEX, hispanic = HISP, edu = SCHL,
           homeown = TEN, st_fips = ST) %>%
    mutate(across(c(age, hispanic, sex, race_White, race_Black, race_Asian,
                    income, edu, homeown, st_fips, PWGTP), as.numeric)) %>%
    filter(age >= AGE_MIN) %>%
    mutate(
      eth  = ifelse(hispanic != 1, "Hispanic", "Non-Hispanic"),
      race = case_when(race_White == 1 ~ "White", race_Black == 1 ~ "Black",
                       race_Asian == 1 ~ "Asian", TRUE ~ "Others"),
      sex  = case_when(sex == 1 ~ "Male", sex == 2 ~ "Female", TRUE ~ NA_character_),
      race_eth = case_when(
        eth == "Hispanic"                       ~ "Hispanic",
        race == "White" & eth == "Non-Hispanic" ~ "NH White",
        race == "Black" & eth == "Non-Hispanic" ~ "NH Black",
        race == "Asian" & eth == "Non-Hispanic" ~ "NH Asian",
        TRUE                                    ~ "Others"),
      income = case_when(
        income >= -60000 & income < 35000  ~ "<35k",
        income >= 35000  & income < 50000  ~ "35k-50k",
        income >= 50000  & income < 75000  ~ "50k-75k",
        income >= 75000  & income < 100000 ~ "75k-100k",
        income >= 100000                   ~ ">100k",
        TRUE                               ~ NA_character_),
      edu = case_when(
        edu < 12              ~ "Less than highschool",
        edu >= 12 & edu < 16  ~ "Some highschool",
        edu == 16 | edu == 17 ~ "Highschool graduate",
        edu >= 18 & edu < 21  ~ "Some college",
        edu >= 21             ~ "College graduate or advanced",
        TRUE                  ~ NA_character_),
      homeown = case_when(
        homeown == 1 | homeown == 2 ~ "Own",
        homeown == 3                ~ "Rent",
        homeown == 4                ~ "Others",
        TRUE                        ~ NA_character_),
      state    = unname(fips2usps[sprintf("%02d", st_fips)]),
      region   = unname(state_region_map[state]),      # official Census map, as on the AoU side
      agegroup = make_agegroup(age)
    ) %>%
    select(state, region, PWGTP, age, agegroup, sex, race_eth, income, edu, homeown)

  message("PUMS ", AGE_MIN, "+ records: ", nrow(pums40), " | 40+ population (sum PWGTP): ",
          format(round(sum(pums40$PWGTP)), big.mark = ","))
  if (any(is.na(pums40$state))) warning("Unmapped ST codes in PUMS — check fips2usps.")
  write.csv(pums40, F_PUMS, row.names = FALSE)
  message("Wrote ", F_PUMS)
} else {
  message("Stage 1 cached: ", F_PUMS)
}

########################################################################
## STAGE 2 — AoU 40+ cohort + Cataract indicator (Workbench 2.0 / BigQuery)
## Same CDR-resolution grammar, cohort (has_ehr_data = 1), survey concept
## ids, and recodes as 02_AoU_Prep.R / Build_Disease_Indicators.R. Adds the
## phenotype query: number of DISTINCT dates carrying a matching ICD source
## code, so the case threshold (MIN_CODE_DATES) can be changed in stage 4
## without re-querying.
## No setup is needed in Researcher Workbench 2.0 when the workspace has the
## "All of Us Controlled Tier" data collection attached: the CDR is found
## through the `wb` CLI (see "Locating the CDR" below). Optional overrides:
##   Sys.setenv(AOU_CDR_RESOURCE_ID = "<resource id from wb resource list>")
##   Sys.setenv(AOU_CDR_DATASET     = "project.dataset")
########################################################################

if (!stage_cached(F_AOU)) {

  suppressPackageStartupMessages({ library(bigrquery); library(lubridate) })

  ## ------------------------------------------------------------------
  ## Locating the CDR in Researcher Workbench 2.0 (Verily Workbench)
  ##
  ## The CDR is no longer announced through WORKSPACE_CDR. It is a BigQuery
  ## dataset RESOURCE attached to the workspace (the "All of Us Controlled
  ## Tier" data collection), and the Workbench CLI `wb` — installed and
  ## authenticated in every cloud app — is the way to find it:
  ##   wb resource list --type=BQ_DATASET --format=JSON   -> the workspace's datasets
  ##   wb resolve --id=<resource id>                      -> "project.dataset"
  ## (The WORKBENCH_<resource> variables the docs mention are only set inside
  ## `wb` subprocesses, so Sys.getenv() cannot see them from R.) The billing
  ## project is the workspace's Google project, GOOGLE_CLOUD_PROJECT.
  ##
  ## Resolution order: AOU_CDR_DATASET (explicit "project.dataset" override)
  ## > AOU_CDR_RESOURCE_ID (explicit resource id) > automatic discovery from
  ## `wb resource list` > legacy WORKSPACE_CDR. Nothing needs to be set when
  ## the workspace has exactly one CDR data collection attached.
  ## ------------------------------------------------------------------
  library(jsonlite)

  AOU_CDR_RESOURCE_ID <- Sys.getenv("AOU_CDR_RESOURCE_ID", unset = "")
  AOU_CDR_DATASET     <- Sys.getenv("AOU_CDR_DATASET",     unset = "")
  BQ_BILLING_PROJECT  <- Sys.getenv("BQ_BILLING_PROJECT",  unset = "")

  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
  is_nonempty <- function(x) length(x) == 1 && !is.na(x) && nzchar(trimws(x))
  run_cmd <- function(cmd, args) {
    out <- tryCatch(suppressWarnings(system2(cmd, args = args, stdout = TRUE, stderr = TRUE)),
                    error = function(e) character())
    out[!grepl("^\\s*$", out)]
  }
  ## Parse the JSON a `wb ... --format=JSON` call prints (skipping any notice
  ## lines the CLI writes before the JSON body)
  wb_json <- function(args) {
    out <- run_cmd("wb", args)
    if (!length(out)) return(NULL)
    txt   <- paste(out, collapse = "\n")
    start <- regexpr("[\\[{]", txt)
    if (start < 0) return(NULL)
    tryCatch(fromJSON(substr(txt, start, nchar(txt)), simplifyVector = FALSE),
             error = function(e) NULL)
  }
  resolve_wb_resource <- function(resource_id) {
    for (args in list(c("resolve", paste0("--id=", resource_id), "--bq-path=FULL_PATH"),
                      c("resolve", paste0("--id=", resource_id)),
                      c("resource", "resolve", paste0("--id=", resource_id)))) {
      out <- run_cmd("wb", args)
      hit <- out[grepl("^[A-Za-z0-9_.-]+\\.[A-Za-z0-9_]+$", trimws(out))]
      if (length(hit)) return(trimws(hit[1]))
    }
    stop("Could not resolve Workbench resource '", resource_id, "' with `wb resolve`.")
  }
  ## Automatic discovery: the BigQuery dataset resources of this workspace
  list_bq_resources <- function() {
    res <- wb_json(c("resource", "list", "--type=BQ_DATASET", "--format=JSON"))
    if (is.null(res)) return(NULL)
    if (!is.null(names(res)) && !is.null(res$id)) res <- list(res)     # single object
    tab <- bind_rows(lapply(res, function(r) tibble(
      id          = as.character(r$id %||% r$name %||% NA_character_),
      stewardship = as.character(r$stewardshipType %||% r$stewardship %||% NA_character_),
      description = as.character(r$description %||% ""),
      projectId   = as.character(r$projectId %||% r$resourceAttributes$gcpBqDataset$projectId %||% NA_character_),
      datasetId   = as.character(r$datasetId %||% r$resourceAttributes$gcpBqDataset$datasetId %||% NA_character_)
    )))
    if (nrow(tab) == 0) NULL else tab
  }
  get_cdr_dataset <- function() {
    if (is_nonempty(AOU_CDR_DATASET))     return(AOU_CDR_DATASET)
    if (is_nonempty(AOU_CDR_RESOURCE_ID)) return(resolve_wb_resource(AOU_CDR_RESOURCE_ID))

    tab <- list_bq_resources()
    if (!is.null(tab)) {
      message("BigQuery dataset resources in this workspace:\n",
              paste(sprintf("  %-45s %-11s %s", tab$id, coalesce(tab$stewardship, "?"),
                            ifelse(is.na(tab$datasetId), "", paste0(tab$projectId, ".", tab$datasetId))),
                    collapse = "\n"))
      ## The CDR is a REFERENCED data collection; datasets the workspace itself
      ## created (scratch results, Data Explorer exports) are CONTROLLED.
      cand <- tab %>% filter(is.na(stewardship) | toupper(stewardship) != "CONTROLLED")
      if (nrow(cand) == 0) cand <- tab
      looks_cdr <- grepl("cdr|controlled|registered|all.?of.?us|aou|[CR]20[0-9]{2}Q[1-4]R[0-9]+",
                         paste(cand$id, cand$description, cand$datasetId), ignore.case = TRUE)
      if (any(looks_cdr) && sum(looks_cdr) < nrow(cand)) cand <- cand[looks_cdr, ]
      if (nrow(cand) > 1) {
        ct <- grepl("controlled|_ct|ct_", paste(cand$id, cand$description, cand$datasetId), ignore.case = TRUE)
        if (sum(ct) == 1) cand <- cand[ct, ]
      }
      if (nrow(cand) == 1) {
        message("Using data collection resource '", cand$id, "'")
        if (!is.na(cand$projectId) && !is.na(cand$datasetId))
          return(paste0(cand$projectId, ".", cand$datasetId))
        return(resolve_wb_resource(cand$id))
      }
      stop("Several BigQuery dataset resources could be the CDR: ",
           paste(cand$id, collapse = ", "),
           ".\nPick one with  Sys.setenv(AOU_CDR_RESOURCE_ID = '<resource id>')  and re-run.")
    }

    legacy <- Sys.getenv("WORKSPACE_CDR", unset = "")
    if (is_nonempty(legacy)) return(legacy)
    stop("No AoU CDR dataset found: `wb resource list` returned no BigQuery dataset ",
         "resources and WORKSPACE_CDR is unset.\n",
         "Attach the 'All of Us Controlled Tier' data collection to this workspace ",
         "(Resources > Add data collection), or set\n",
         "  Sys.setenv(AOU_CDR_DATASET = 'project.dataset')   # e.g. the id shown by `wb resource describe`")
  }
  parse_bq_dataset <- function(dataset_path) {
    parts <- strsplit(gsub("`", "", trimws(dataset_path)), "\\.")[[1]]
    if (length(parts) != 2) stop("Expected 'project.dataset', got: ", dataset_path)
    list(project = parts[1], dataset = parts[2], path = paste(parts, collapse = "."))
  }
  get_billing_project <- function(default_project) {
    cand <- c(BQ_BILLING_PROJECT,
              Sys.getenv("GOOGLE_CLOUD_PROJECT", unset = ""),   # Workbench 2.0 workspace project
              Sys.getenv("GOOGLE_PROJECT",       unset = ""),   # legacy Workbench
              Sys.getenv("GCLOUD_PROJECT",       unset = ""))
    cand <- cand[nzchar(cand)]
    if (length(cand) > 0) return(cand[1])
    ws <- wb_json(c("workspace", "describe", "--format=JSON"))
    wp <- ws$googleProjectId %||% ws$gcpProjectId %||% ws$projectId %||% NULL
    if (is_nonempty(wp)) return(wp)
    gp <- run_cmd("gcloud", c("config", "get-value", "project"))
    gp <- gp[nzchar(gp) & !grepl("^Your active configuration", gp)]
    if (length(gp) > 0 && !grepl("ERROR|unset", gp[1], ignore.case = TRUE)) return(gp[1])
    default_project
  }

  cdr             <- parse_bq_dataset(get_cdr_dataset())
  billing_project <- get_billing_project(cdr$project)
  message("Using AoU CDR dataset: ", cdr$path, " | billing project: ", billing_project)

  fq_table <- function(table_name) paste0("`", cdr$path, ".", table_name, "`")
  table_exists <- function(table_name) {
    isTRUE(tryCatch(bq_table_exists(bq_table(cdr$project, cdr$dataset, table_name)),
                    error = function(e) FALSE))
  }
  ## person_id as double throughout (joins against read_csv output downstream)
  run_bq <- function(sql, bigint = "numeric") {
    job <- bq_project_query(x = billing_project, query = sql, use_legacy_sql = FALSE)
    bq_table_download(job, bigint = bigint, quiet = TRUE)
  }

  ## --- Cohort: participants with EHR data (as in 02_AoU_Prep.R) ---
  if (table_exists("cb_search_person")) {
    ehr_cohort_sql <- paste0("SELECT DISTINCT person_id FROM ", fq_table("cb_search_person"),
                             " WHERE has_ehr_data = 1")
  } else {
    warning("cb_search_person not found — using ALL persons (not the EHR-only cohort).")
    ehr_cohort_sql <- paste0("SELECT DISTINCT person_id FROM ", fq_table("person"))
  }

  dataset_person_sql <- paste0("
    SELECT person.person_id,
           p_gender_concept.concept_name       AS gender,
           person.birth_datetime               AS date_of_birth,
           p_race_concept.concept_name         AS race,
           p_ethnicity_concept.concept_name    AS ethnicity,
           p_sex_at_birth_concept.concept_name AS sex_at_birth
    FROM ", fq_table("person"), " person
    LEFT JOIN ", fq_table("concept"), " p_gender_concept
      ON person.gender_concept_id = p_gender_concept.concept_id
    LEFT JOIN ", fq_table("concept"), " p_race_concept
      ON person.race_concept_id = p_race_concept.concept_id
    LEFT JOIN ", fq_table("concept"), " p_ethnicity_concept
      ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id
    LEFT JOIN ", fq_table("concept"), " p_sex_at_birth_concept
      ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
    WHERE person.person_id IN (", ehr_cohort_sql, ")")

  dataset_survey_sql <- paste0("
    SELECT answer.person_id, answer.question_concept_id, answer.question, answer.answer
    FROM ", fq_table("ds_survey"), " answer
    WHERE question_concept_id IN (1585370, 1585375, 1585899, 1585940, 43530593)
      AND answer.person_id IN (", ehr_cohort_sql, ")")

  dataset_state_sql <- paste0("
    SELECT person.person_id, person.state_of_residence_source_value AS state
    FROM ", fq_table("person_ext"), " person
    WHERE person.person_id IN (", ehr_cohort_sql, ")")

  ## --- Phenotype: distinct dates with a matching ICD-9/10 SOURCE code ---
  ## Matched on the source concept's concept_code (dotted ICD code) OR the raw
  ## source value, by prefix, so 7-character ICD-10 codes are captured and
  ## neighbouring codes are not. Both condition_occurrence and observation
  ## are searched (OMOP routes ICD "status" codes such as Z96.1 / V43.1 to the
  ## observation table, and some sites record diagnoses there); the dates are
  ## pooled across the two tables.
  pheno_regex <- paste0("^(", paste(gsub(".", "[.]", c(PHENO_ICD9, PHENO_ICD10), fixed = TRUE),
                                    collapse = "|"), ")")
  has_obs <- table_exists("observation")
  if (!has_obs) warning("observation table not found — status codes such as Z96.1 cannot be counted.")
  src_sql <- function(tab, scid, sval, dcol) paste0("
    SELECT t.person_id, t.", dcol, " AS d
    FROM ", fq_table(tab), " t
    JOIN ", fq_table("concept"), " mc ON t.", scid, " = mc.concept_id
    WHERE mc.vocabulary_id IN ('ICD9CM', 'ICD10CM')
      AND (", code_cond("mc.concept_code", pheno_regex), "
           OR ", code_cond(paste0("IFNULL(t.", sval, ", '')"), pheno_regex), ")
      AND t.person_id IN (", ehr_cohort_sql, ")")
  parts <- src_sql("condition_occurrence", "condition_source_concept_id", "condition_source_value", "condition_start_date")
  if (has_obs)
    parts <- paste(parts, "UNION ALL",
                   src_sql("observation", "observation_source_concept_id", "observation_source_value", "observation_date"))
  pheno_sql <- paste0("SELECT person_id, COUNT(DISTINCT d) AS n_", PHENO_ID, "_dates FROM (", parts, ") GROUP BY person_id")
  message(PHENO_LABEL, " ICD regex: ", pheno_regex)
  message(PHENO_LABEL, " case definition: ", PHENO_DEF,
          if (PHENO_DEF == "phecodeX") paste0(" (PhecodeX ", PHECODEX_CATARACT$phecode, ", ", length(PHECODEX_CODES),
                                               " exact ICD codes)") else " (VEHSS crosswalk prefixes)")

  person_df <- run_bq(dataset_person_sql) %>%
    mutate(person_id = as.numeric(person_id),
           across(c(gender, race, ethnicity, sex_at_birth), as.character))
  survey_df <- run_bq(dataset_survey_sql) %>%
    mutate(person_id = as.numeric(person_id), across(c(question, answer), as.character))
  state_df  <- run_bq(dataset_state_sql) %>%
    mutate(person_id = as.numeric(person_id), state = as.character(state))
  pheno_df  <- run_bq(pheno_sql) %>%
    mutate(person_id = as.numeric(person_id), n_cataract_dates = as.numeric(n_cataract_dates))
  message("Rows: person=", nrow(person_df), ", survey=", nrow(survey_df),
          ", state=", nrow(state_df), ", persons with any ", PHENO_LABEL, " code=", nrow(pheno_df))

  ## --- Recode (identical to 02_AoU_Prep.R), then restrict to age >= AGE_MIN ---
  survey_wide <- survey_df %>%
    select(-question_concept_id) %>%
    filter(question != "The Basics: Sexual Orientation") %>%
    pivot_wider(names_from = question, values_from = answer, values_fn = first)

  dt0 <- person_df %>%
    full_join(survey_wide, by = "person_id") %>%
    full_join(state_df,    by = "person_id") %>%
    rename(income  = `Income: Annual Income`,
           homeown = `Home Own: Current Home Own`,
           edu     = `Education Level: Highest Grade`) %>%
    mutate(date_of_birth = as.Date(date_of_birth),
           age = round(interval(date_of_birth, AGE_REFERENCE_DATE) /
                         duration(num = 1, units = "years"), 0)) %>%
    mutate(race = case_when(
      race %in% c("I prefer not to answer", "PMI: Skip") ~ NA_character_,
      race %in% c("None Indicated", "None of these", "Middle Eastern or North African",
                  "Native Hawaiian or Other Pacific Islander",
                  "More than one population")             ~ "Others",
      TRUE                                                ~ race)) %>%
    mutate(ethnicity = case_when(
      ethnicity == "Hispanic or Latino"                                     ~ "Yes",
      ethnicity %in% c("No matching concept", "Not Hispanic or Latino",
                       "What Race Ethnicity: Race Ethnicity None Of These") ~ "No",
      ethnicity %in% c("PMI: Prefer Not To Answer", "PMI: Skip")            ~ NA_character_,
      TRUE                                                                  ~ ethnicity)) %>%
    mutate(sex = case_when(
      sex_at_birth %in% c("I prefer not to answer", "Intersex", "No matching concept",
                          "None", "PMI: Skip") ~ NA_character_,
      TRUE                                     ~ sex_at_birth)) %>%
    mutate(income = case_when(
      income %in% c("Annual Income: 100k 150k", "Annual Income: 150k 200k",
                    "Annual Income: more 200k")               ~ ">100k",
      income == "Annual Income: 75k 100k"                     ~ "75k-100k",
      income == "Annual Income: 50k 75k"                      ~ "50k-75k",
      income == "Annual Income: 35k 50k"                      ~ "35k-50k",
      income %in% c("Annual Income: 25k 35k", "Annual Income: 10k 25k",
                    "Annual Income: less 10k")                ~ "<35k",
      TRUE                                                    ~ NA_character_)) %>%
    mutate(homeown = case_when(
      homeown %in% c("Current Home Own: Other Arrangement", "PMI: Dont Know") ~ "Others",
      homeown == "Current Home Own: Own"                                      ~ "Own",
      homeown == "Current Home Own: Rent"                                     ~ "Rent",
      TRUE                                                                    ~ NA_character_)) %>%
    mutate(edu = case_when(
      edu %in% c("Highest Grade: Advanced Degree",
                 "Highest Grade: College Graduate")   ~ "College graduate or advanced",
      edu == "Highest Grade: College One to Three"    ~ "Some college",
      edu == "Highest Grade: Twelve Or GED"           ~ "Highschool graduate",
      edu == "Highest Grade: Nine Through Eleven"     ~ "Some highschool",
      edu %in% c("Highest Grade: Never Attended", "Highest Grade: One Through Four",
                 "Highest Grade: Five Through Eight") ~ "Less than highschool",
      TRUE                                            ~ NA_character_)) %>%
    mutate(state = ifelse(str_detect(state, "^PII State: "),
                          str_extract(state, "(?<=PII State: )\\w{2}"), NA_character_))

  aou40 <- dt0 %>%
    filter(!is.na(age), age >= AGE_MIN) %>%
    mutate(race_eth = case_when(
      ethnicity == "Yes"                       ~ "Hispanic",
      race      == "Asian"                     ~ "NH Asian",
      race      == "Black or African American" ~ "NH Black",
      race      == "White"                     ~ "NH White",
      is.na(race) | is.na(ethnicity)           ~ NA_character_,
      TRUE                                     ~ "Others"),
      agegroup = make_agegroup(age),
      region   = unname(state_region_map[state])) %>%
    select(person_id, state, region, age, agegroup, sex, race_eth, income, edu, homeown) %>%
    harmonize_factors() %>%                      # drops rows missing any covariate / region
    left_join(pheno_df, by = "person_id") %>%
    mutate(n_cataract_dates = coalesce(n_cataract_dates, 0))

  message("AoU ", AGE_MIN, "+ participants with complete covariates + state: ", nrow(aou40),
          " | with >= ", MIN_CODE_DATES, " ", PHENO_LABEL, " code dates: ", sum(aou40$n_cataract_dates >= MIN_CODE_DATES))
  write_excel_csv(aou40, F_AOU)
  message("Wrote ", F_AOU)
} else {
  message("Stage 2 cached: ", F_AOU)
}

########################################################################
## STAGE 3 — Global RAILS + Subgroup RAILS (region) on the 40+ cell tables
## Aggregation as in 01_PUMS_Prep.R / 02_AoU_Prep.R; estimation as in
## 04_Global_RAILS.R (7 covariates incl. region) and 09_Sub_RAILS_region.R
## (6 covariates within each region, region's own PUMS total as nsiz).
########################################################################

if (!stage_cached(F_WEIGHTS)) {

  source(SUB_FUN_FILE, chdir = TRUE)   # fun.sub.rails.threeway (+ AoU_Fun.R next to it)

  dt_aou <- read_csv(F_AOU, show_col_types = FALSE) %>%
    mutate(person_id = as.numeric(person_id)) %>%
    harmonize_factors() %>%
    mutate(weight = 1)

  ## ---- Reference weights: restrict to 40+ (done in stage 1), then RESCALE ----
  ## harmonize_factors() drops PUMS records with a missing covariate (income /
  ## edu / homeown can be NA) and records outside state_region_map (none: DC is mapped to the South). The
  ## surviving weights are rescaled so they sum to the 40+ target population,
  ## NOT left at their raw (reduced) sum and NOT at the all-adult total the
  ## main pipeline uses. The anchor is taken over every 40+ record of the
  ## mapped states before that drop, mirroring n_adult in 01_PUMS_Prep.R.
  pums40 <- read_csv(F_PUMS, show_col_types = FALSE) %>%
    filter(age >= AGE_MIN)                              # no-op unless F_PUMS predates AGE_MIN
  n_pop40_all    <- sum(pums40$PWGTP, na.rm = TRUE)                       # all 40+ incl. DC
  n_pop40_mapped <- sum(pums40$PWGTP[!is.na(pums40$region)], na.rm = TRUE) # 40+ in mapped states
  n_pop40 <- if (is.finite(TARGET_POP)) TARGET_POP else n_pop40_mapped
  dt_pums <- pums40 %>%
    harmonize_factors() %>%
    mutate(weight = PWGTP / sum(PWGTP) * n_pop40)      # rescale to the 40+ target
  message("PUMS ", AGE_MIN, "+ population: all states+DC = ", format(round(n_pop40_all), big.mark = ","),
          " | mapped states = ", format(round(n_pop40_mapped), big.mark = ","),
          " | complete-covariate raw sum = ",
          format(round(sum(dt_pums$PWGTP)), big.mark = ","),
          "\n  -> target population (nsiz) = ", format(round(n_pop40), big.mark = ","),
          if (is.finite(TARGET_POP)) "  [TARGET_POP override]" else "  [ACS 40+, mapped states]",
          "\n  VEHSS-implied 2019 40+ population (US row) = ",
          format(round(vehss_us$pop_vehss), big.mark = ","))
  stopifnot(abs(sum(dt_pums$weight) - n_pop40) < 1e-6 * n_pop40)

  cat_formula <- formula(paste0("weight ~ ", paste(NAMES_7, collapse = "+")))
  dt_agg_aou  <- aggregate(cat_formula, data = dt_aou,  sum)
  dt_agg_pums <- aggregate(cat_formula, data = dt_pums, sum)
  message("Cells: AoU=", nrow(dt_agg_aou), " (n=", sum(dt_agg_aou$weight), "), PUMS=",
          nrow(dt_agg_pums), " (N=", format(round(sum(dt_agg_pums$weight)), big.mark = ","), ")")
  print(dt_agg_pums %>% group_by(region) %>%
          summarise(pop_40plus = round(sum(weight)), .groups = "drop") %>%
          mutate(share = round(pop_40plus / sum(pop_40plus), 4)))

  ## Population margins: one-, two-, three-way over the 7 covariates
  twovars     <- combn(NAMES_7, 2, FUN = function(x) paste(x, collapse = ":"))
  threevars   <- combn(NAMES_7, 3, FUN = function(x) paste(x, collapse = ":"))
  max_formula <- formula(paste0("~", paste(c(NAMES_7, twovars, threevars), collapse = "+")))
  mat_max     <- sparse.model.matrix(max_formula, data = dt_agg_pums, keep.order = TRUE)
  pop_totals  <- setNames(as.numeric(Matrix::crossprod(mat_max, dt_agg_pums$weight)),
                          colnames(mat_max))
  if (any(pop_totals == 0))
    warning(sum(pop_totals == 0), " zero PUMS margins — raking on terms involving them will fail.")

  ## Global RAILS (region is a calibration variable)
  message("=== Global RAILS, ", AGE_MIN, "+ ===")
  result_global <- fun.rails.threeway(
    dt_agg_aou = dt_agg_aou, dt_agg_pums = dt_agg_pums, pop_totals = pop_totals,
    names_univar = NAMES_7, alpha = ALPHA
  )
  message("Global calibrated model: ", result_global$calibrated_terms[1])

  ## Subgroup RAILS by region (sex stays in the model, region is the stratum).
  ## Within a region the two-way base model of fun.nps can hit a singular
  ## Newton step when the region's AoU or PUMS cells are too sparse (an empty
  ## one-/two-way cell is the usual cause, but near-separation of a sparse cell
  ## pattern can do it too). Report empty cells first, then run each region
  ## under tryCatch so one failing region yields NA w_subrails (and a warning)
  ## instead of aborting the whole run. If a region fails, coarsen AGE_BREAKS
  ## (e.g. c(40, 65, Inf)) and re-run (agegroup is re-cut here in stage 3).
  message("=== Subgroup RAILS by region, ", AGE_MIN, "+ ===")
  zero_margins <- function(agg, vars) {
    f <- formula(paste0("~", paste(c(vars, combn(vars, 2, FUN = function(x) paste(x, collapse = ":"))),
                                   collapse = "+")))
    m <- sparse.model.matrix(f, agg)
    colnames(m)[as.numeric(Matrix::crossprod(m, agg$weight)) == 0]
  }
  for (lev in REGION_LEVELS) {
    z_aou  <- zero_margins(dt_agg_aou  %>% filter(region == lev), NAMES_6)
    z_pums <- zero_margins(dt_agg_pums %>% filter(region == lev), NAMES_6)
    if (length(z_aou) || length(z_pums))
      warning("Region ", lev, ": empty one-/two-way cells — AoU: ",
              paste(z_aou, collapse = ", "), " | PUMS: ", paste(z_pums, collapse = ", "),
              " (the within-region two-way model may not be estimable)", immediate. = TRUE)
  }
  result_region <- bind_rows(lapply(REGION_LEVELS, function(lev) {
    tryCatch(
      fun.sub.rails.threeway(
        dt_agg_aou = dt_agg_aou %>% filter(region == lev),   # factor levels kept
        dt_agg_pums = dt_agg_pums %>% filter(region == lev),
        subgroup_var = "region", names_univar = NAMES_6, alpha = ALPHA),
      error = function(e) {
        warning("Subgroup RAILS FAILED for region ", lev, ": ", conditionMessage(e),
                " — w_subrails will be NA for that region.", immediate. = TRUE)
        NULL
      })
  }))
  if (nrow(result_region) == 0) stop("Subgroup RAILS failed in every region.")
  print(result_region %>% group_by(subgroup_run) %>%
          summarise(n_cells = n(), total_srails = sum(d_rails * weight, na.rm = TRUE),
                    calibrated = first(calibrated_terms), .groups = "drop"))

  write_excel_csv(result_global, F_CELLS_G)
  write_excel_csv(result_region, F_CELLS_S)

  ## Per-individual weights joined back (d_* are already per-individual)
  dt_w <- dt_aou %>%
    left_join(
      result_global %>%
        transmute(across(all_of(NAMES_7)),
                  w_unweighted = d_unweighted, w_cal1 = d_cal1, w_cal2 = d_cal2,
                  w_nps1 = d_nps1, w_nps2 = d_nps2, w_nps1_rake = d_nps1_rake,
                  w_nps2_rake = d_nps2_rake, w_rails = d_rails,
                  global_calibrated_terms = calibrated_terms),
      by = NAMES_7) %>%
    left_join(
      result_region %>%
        transmute(across(all_of(NAMES_7)), w_subrails = d_rails,
                  region_calibrated_terms = calibrated_terms),
      by = NAMES_7) %>%
    select(-weight)

  if (all(is.na(dt_w$w_rails)))    stop("Global RAILS produced no weights (no converged three-way step).")
  if (all(is.na(dt_w$w_subrails))) stop("Subgroup RAILS produced no weights in any region.")
  ## FALLBACK: a region whose subgroup fit did not converge keeps the Global
  ## RAILS weight as its S-RAILS weight (flagged in subrails_source), so the
  ## S-RAILS estimates stay complete; the user is told loudly here and again
  ## in stage 4, and the affected regions are named on the figures.
  miss_sub <- as.character(dt_w %>% filter(is.na(w_subrails)) %>% distinct(region) %>% pull(region))
  dt_w <- dt_w %>%
    mutate(subrails_source = ifelse(is.na(w_subrails), "global_fallback", "subgroup"),
           w_subrails      = coalesce(w_subrails, w_rails))
  if (length(miss_sub)) {
    msg <- paste0("SUBGROUP RAILS DID NOT CONVERGE IN: ", paste(miss_sub, collapse = ", "),
                  ". Global RAILS weights are used as the S-RAILS weights there ",
                  "(subrails_source = 'global_fallback').")
    message("\n!!! ", msg, "\n"); warning(msg, immediate. = TRUE)
  }

  ## Every weight column must sum to the 40+ target (S-RAILS: region by region,
  ## after the fallback, so a fallback region shows its Global RAILS total)
  message("Weight totals: G-RAILS=", format(round(sum(dt_w$w_rails, na.rm = TRUE)), big.mark = ","),
          " | S-RAILS=", format(round(sum(dt_w$w_subrails, na.rm = TRUE)), big.mark = ","),
          " | target 40+ population=", format(round(n_pop40), big.mark = ","))
  if (abs(sum(dt_w$w_rails, na.rm = TRUE) - n_pop40) > 1e-3 * n_pop40)
    warning("G-RAILS weights do not sum to the 40+ target population — check the calibration.")
  print(dt_w %>% group_by(region) %>%
          summarise(sum_w_rails = round(sum(w_rails, na.rm = TRUE)),
                    sum_w_subrails = round(sum(w_subrails, na.rm = TRUE)), .groups = "drop") %>%
          left_join(dt_agg_pums %>% group_by(region) %>%
                      summarise(pums_40plus = round(sum(weight)), .groups = "drop"), by = "region"))
  write_excel_csv(dt_w, F_WEIGHTS)
  message("Wrote ", F_WEIGHTS)
} else {
  message("Stage 3 cached: ", F_WEIGHTS)
}

########################################################################
## STAGE 4 — Cataract prevalence vs VEHSS: overall and by Age x Sex x Race
########################################################################

df <- read_csv(F_WEIGHTS, show_col_types = FALSE) %>%
  mutate(y = as.integer(n_cataract_dates >= MIN_CODE_DATES), one = 1) %>%
  filter(!is.na(state), !is.na(region), !is.na(w_rails))
if (!"subrails_source" %in% names(df))          # weights file from an older run
  df <- df %>% mutate(subrails_source = ifelse(is.na(w_subrails), "global_fallback", "subgroup"),
                      w_subrails = coalesce(w_subrails, w_rails))
fallback_regions <- sort(unique(as.character(df$region[df$subrails_source == "global_fallback"])))
fallback_note <- if (length(fallback_regions))
  paste0("S-RAILS = Global RAILS fallback in: ", paste(fallback_regions, collapse = ", ")) else NULL
REF_NOTE      <- "Medicare FFS + Advantage claims 2022: diagnosed cataract, any VEHSS cataract subgroup incl. pseudophakia; denominator = Medicare beneficiaries (40-64 = disability / ESRD beneficiaries only), not the resident population"
fallback_note <- c(paste0("Reference: ", REF_NOTE), fallback_note)   # shown on every figure
if (length(fallback_regions))
  message("\n!!! NOTE: Subgroup RAILS did not converge in ", paste(fallback_regions, collapse = ", "),
          " — Global RAILS weights stand in for S-RAILS there (S-RAILS == G-RAILS in those states).\n")
message("Participants ", AGE_MIN, "+ with weights: ", nrow(df),
        " | ", PHENO_LABEL, " cases: ", sum(df$y), " (", sprintf("%.2f%%", 100 * mean(df$y)), " unweighted)")

## ---- VEHSS-aligned domain variables (AoU and PUMS) ----------------------
## VEHSS race groups: Black / Hispanic / White (non-Hispanic) / Other; the
## pipeline's NH Asian and Others both fall into VEHSS "Other".
to_race_v <- function(race_eth) case_when(
  race_eth == "Hispanic" ~ "Hispanic", race_eth == "NH Black" ~ "Black",
  race_eth == "NH White" ~ "White",    TRUE                   ~ "Other")
df <- df %>% mutate(age_v = as.character(agegroup), sex_v = as.character(sex),
                    race_v = to_race_v(race_eth))

pums40 <- read_csv(F_PUMS, show_col_types = FALSE) %>%
  filter(age >= AGE_MIN, !is.na(state)) %>%
  mutate(age_v = make_agegroup(age), sex_v = as.character(sex), race_v = to_race_v(race_eth))

RUN_COMBOS <- identical(AGE_LABELS, c("40-64", "65-84", "85+"))
if (!RUN_COMBOS)
  warning("AGE_BREAKS are not the VEHSS strata (40-64 / 65-84 / 85+): the age-specific ",
          "combinations are skipped; only the age = 40+ combinations are compared.", immediate. = TRUE)

## Subset to one Age x Sex x Race selection ("40+" / "Both" / "All" = no restriction)
dom <- function(d, A, S, R) {
  d %>% filter((A == "40+" | age_v == A), (S == "Both" | sex_v == S), (R == "All" | race_v == R))
}

METHODS <- c(unweighted = "one", grails = "w_rails", srails = "w_subrails")
METHOD_LABELS <- c(vehss = "VEHSS Medicare FFS+MA diagnosed (2022)", unweighted = "AoU unweighted",
                   grails = "AoU G-RAILS", srails = "AoU S-RAILS (region)")
METHOD_COLS   <- c(vehss = "black", unweighted = "#E08214", grails = "#2166AC", srails = "#B2182B")

## Weighted prevalence with a linearized 95% CI (as 05_Prevalence_Analysis.R)
wt.prev <- function(y, w) {
  ok <- !is.na(y) & !is.na(w); y <- y[ok]; w <- w[ok]
  if (length(y) == 0 || sum(w) <= 0) return(c(est = NA_real_, lb = NA_real_, ub = NA_real_))
  p <- sum(w * y) / sum(w)
  v <- sum(w^2 * (y - p)^2) / sum(w)^2
  c(est = p, lb = p - qnorm(0.975) * sqrt(v), ub = p + qnorm(0.975) * sqrt(v))
}

## Long table: one row per group x method (n, n_cases, est, lb, ub)
prev_table <- function(d, by = character(0)) {
  if (nrow(d) == 0) return(tibble(n = integer(), n_cases = integer(), est = numeric(),
                                  lb = numeric(), ub = numeric(), method = factor(character(), levels = names(METHODS))))
  bind_rows(lapply(names(METHODS), function(m) {
    d %>%
      mutate(w = .data[[METHODS[m]]]) %>%
      group_by(across(all_of(by))) %>%
      summarise(n = sum(!is.na(w)), n_cases = sum(y[!is.na(w)]),
                est = wt.prev(y, w)[["est"]],
                lb  = wt.prev(y, w)[["lb"]],
                ub  = wt.prev(y, w)[["ub"]],
                .groups = "drop") %>%
      mutate(method = m)
  })) %>%
    mutate(method = factor(method, levels = names(METHODS)))
}

## Manuscript application metrics with z = state (10_State_Prevalence_Maps.R)
H_bern  <- function(p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); -p * log(p) - (1 - p) * log(1 - p) }
js_bern <- function(pg, ps) H_bern((pg + ps) / 2) - 0.5 * H_bern(pg) - 0.5 * H_bern(ps)
## rho is renormalized over the states that enter (un-suppressed), so that
## combinations with many suppressed states are not pushed towards D = 0.
weighted_js  <- function(rho, d) { ok <- is.finite(d) & is.finite(rho); if (!any(ok)) return(NA_real_); sum(rho[ok] * d[ok]) / sum(rho[ok]) }
disc_entropy <- function(rho, d) {
  ok  <- is.finite(d) & is.finite(rho)
  c_a <- rho[ok] * d[ok]
  Dj  <- sum(c_a)
  if (length(c_a) <= 1 || Dj <= 0) return(NA_real_)   # < 2 states, or G-RAILS == S-RAILS everywhere
  q <- c_a / Dj; q <- q[q > 0]
  -sum(q * log(q)) / log(length(c_a))
}

########################################################################
## State-level comparison for ONE Age x Sex x Race selection
## Returns the long state table (all VEHSS states x methods) with VEHSS
## values, ratios, S-vs-G divergence, and the selection's summary row.
########################################################################

vehss_states_all <- sort(unique(vehss_all$state[vehss_all$state != "US"]))

compare_combo <- function(A, S, R) {
  d  <- dom(df, A, S, R)
  vs <- vehss_all %>% filter(age == A, sex == S, race == R)
  v_state <- vs %>% filter(state != "US") %>% select(-age, -sex, -race)
  v_us    <- vs %>% filter(state == "US")

  ## PUMS domain population by state -> rho_a and the domain's population
  pdom <- dom(pums40, A, S, R) %>% group_by(state) %>%
    summarise(pop_pums = sum(PWGTP, na.rm = TRUE), .groups = "drop")

  st <- tidyr::expand_grid(state = vehss_states_all, method = names(METHODS)) %>%
    left_join(prev_table(d, "state") %>% mutate(method = as.character(method)),
              by = c("state", "method")) %>%
    left_join(v_state, by = "state") %>%
    left_join(pdom, by = "state") %>%
    mutate(age = A, sex = S, race = R,
           method   = factor(method, levels = names(METHODS)),
           n        = coalesce(n, 0L), n_cases = coalesce(n_cases, 0L),
           suppress = n_cases < SUPPRESS_MIN,
           across(c(est, lb, ub), ~ ifelse(suppress, NA_real_, .x)),
           ratio_to_vehss = est / prev_vehss,
           log2_ratio     = log2(ratio_to_vehss),
           in_vehss_ci    = !is.na(est) & !is.na(lb_vehss) & est >= lb_vehss & est <= ub_vehss,
           ci_overlap     = !is.na(est) & !is.na(lb_vehss) & lb <= ub_vehss & ub >= lb_vehss)

  ## S-RAILS vs G-RAILS per state: ratio and local JS divergence
  sg <- st %>% select(state, method, est, pop_pums) %>%
    pivot_wider(names_from = method, values_from = est) %>%
    mutate(r_s_g = srails / grails, log2_s_g = log2(r_s_g),
           d_gs  = ifelse(is.na(grails) | is.na(srails), NA_real_, js_bern(grails, srails)),
           rho   = pop_pums / sum(pop_pums, na.rm = TRUE))

  nat <- prev_table(d) %>%
    mutate(prev_vehss = if (nrow(v_us)) v_us$prev_vehss else NA_real_,
           lb_vehss   = if (nrow(v_us)) v_us$lb_vehss   else NA_real_,
           ub_vehss   = if (nrow(v_us)) v_us$ub_vehss   else NA_real_,
           across(c(est, lb, ub), ~ ifelse(n_cases < SUPPRESS_MIN, NA_real_, .x)))

  agree <- function(m) {
    x <- st %>% filter(method == m, !suppress, is.finite(ratio_to_vehss))
    if (nrow(x) < 3) return(tibble(n_states_vehss = nrow(x), gm_ratio = NA_real_,
                                   mean_abs_log_ratio = NA_real_, pct_in_vehss_ci = NA_real_, pearson = NA_real_))
    tibble(n_states_vehss     = nrow(x),
           gm_ratio           = exp(mean(log(x$ratio_to_vehss))),
           mean_abs_log_ratio = mean(abs(log(x$ratio_to_vehss))),
           pct_in_vehss_ci    = 100 * mean(x$in_vehss_ci),
           pearson            = suppressWarnings(cor(x$est, x$prev_vehss)))
  }
  a_s <- agree("srails"); a_g <- agree("grails"); a_u <- agree("unweighted")

  summary <- tibble(
    age = A, sex = S, race = R, combo = paste(A, S, R, sep = " / "),
    is_overall   = (A == "40+" & S == "Both" & R == "All"),
    n_aou        = nrow(d),
    n_cases_aou  = sum(d$y),
    pop_pums     = sum(pdom$pop_pums),
    pop_vehss    = if (nrow(v_us)) v_us$pop_vehss else NA_real_,
    has_vehss    = nrow(v_state) > 0,
    n_states_DE  = sum(is.finite(sg$d_gs)),
    D_js         = weighted_js(sg$rho, sg$d_gs),
    E_js         = disc_entropy(sg$rho, sg$d_gs),
    sum_abs_log_sg = sum(abs(sg$log2_s_g), na.rm = TRUE),
    nat_unweighted = nat$est[nat$method == "unweighted"],
    nat_grails     = nat$est[nat$method == "grails"],
    nat_srails     = nat$est[nat$method == "srails"],
    nat_vehss      = if (nrow(v_us)) v_us$prev_vehss else NA_real_,
    srails_n_states = a_s$n_states_vehss, srails_gm_ratio = a_s$gm_ratio,
    srails_mean_abs_log_ratio = a_s$mean_abs_log_ratio, srails_pct_in_ci = a_s$pct_in_vehss_ci,
    srails_pearson = a_s$pearson,
    grails_gm_ratio = a_g$gm_ratio, grails_mean_abs_log_ratio = a_g$mean_abs_log_ratio,
    grails_pct_in_ci = a_g$pct_in_vehss_ci,
    unw_gm_ratio = a_u$gm_ratio, unw_mean_abs_log_ratio = a_u$mean_abs_log_ratio,
    unw_pct_in_ci = a_u$pct_in_vehss_ci
  )
  list(state = st %>% left_join(sg %>% select(state, r_s_g, log2_s_g, d_gs, rho), by = "state"),
       nat = nat %>% mutate(age = A, sex = S, race = R), summary = summary)
}

########################################################################
## All Age x Sex x Race combinations (the overall one included)
########################################################################

combos <- tidyr::expand_grid(age = AGE_V, sex = SEX_V, race = RACE_V) %>%
  mutate(is_overall = age == "40+" & sex == "Both" & race == "All",
         is_partial = !is_overall & (age == "40+" | sex == "Both" | race == "All"))
if (!RUN_COMBOS) combos <- combos %>% filter(age == "40+")
if (!ALLOW_PARTIAL_COMBOS) combos <- combos %>% filter(!is_partial)
combos <- combos %>% select(age, sex, race)
message("Comparing ", nrow(combos), " Age x Sex x Race combinations (overall + ",
        if (ALLOW_PARTIAL_COMBOS) "all partial and " else "", "fully specified cells) ...")
res <- Map(compare_combo, combos$age, combos$sex, combos$race)

combo_state   <- bind_rows(lapply(res, `[[`, "state"))
combo_nat     <- bind_rows(lapply(res, `[[`, "nat"))
combo_summary <- bind_rows(lapply(res, `[[`, "summary")) %>%
  mutate(age = factor(age, AGE_V), sex = factor(sex, SEX_V), race = factor(race, RACE_V)) %>%
  arrange(age, sex, race)

## ---- Overall (40+ / Both / All): the three legacy tables --------------
ov  <- res[[which(combos$age == "40+" & combos$sex == "Both" & combos$race == "All")]]
state_long <- ov$state %>% mutate(region = unname(state_region_map[state]))
nat <- ov$nat %>% mutate(ratio_to_vehss = est / prev_vehss)
message("\n--- National Cataract prevalence, ", AGE_MIN, "+ (overall) ---")
print(nat %>% transmute(method = METHOD_LABELS[as.character(method)],
                        est = sprintf("%.2f%% (%.2f-%.2f)", 100 * est, 100 * lb, 100 * ub),
                        vehss = sprintf("%.2f%% (%.2f-%.2f)", 100 * prev_vehss, 100 * lb_vehss, 100 * ub_vehss),
                        ratio = round(ratio_to_vehss, 3)))
write_excel_csv(nat, F_NATIONAL)

reg <- prev_table(df, "region") %>%
  left_join(vehss_region, by = "region") %>%
  mutate(ratio_to_vehss = est / prev_vehss, region = factor(region, levels = REGION_LEVELS)) %>%
  arrange(region, method)
message("\n--- Regional Cataract prevalence, ", AGE_MIN, "+ (VEHSS aggregated from states by implied 40+ population) ---")
print(reg %>% transmute(region, method = METHOD_LABELS[as.character(method)], n, n_cases,
                        est = sprintf("%.2f%%", 100 * est), vehss = sprintf("%.2f%%", 100 * prev_vehss),
                        ratio = round(ratio_to_vehss, 3)), n = 20)
write_excel_csv(reg, F_REGION)

message("\nOverall: states suppressed (<", SUPPRESS_MIN, " Cataract cases): ",
        sum(state_long$suppress[state_long$method == "srails"] & state_long$n[state_long$method == "srails"] > 0),
        " | states without AoU rows: ",
        paste(state_long %>% filter(method == "srails", n == 0) %>% pull(state), collapse = ", "))

state_wide <- state_long %>%
  select(state, region, n, n_cases, suppress, method, est, lb, ub, ratio_to_vehss, in_vehss_ci) %>%
  pivot_wider(names_from = method, values_from = c(est, lb, ub, ratio_to_vehss, in_vehss_ci),
              names_glue = "{.value}_{method}") %>%
  left_join(vehss_state %>% select(state, prev_vehss, lb_vehss, ub_vehss, cases_vehss, pop_vehss), by = "state") %>%
  left_join(ov$state %>% distinct(state, r_s_g, log2_s_g, d_gs), by = "state") %>%
  arrange(desc(prev_vehss))
write_excel_csv(state_wide, F_STATE)

summ <- state_long %>%
  filter(!suppress, is.finite(ratio_to_vehss)) %>%
  group_by(method) %>%
  summarise(n_states = n(), mean_ratio = mean(ratio_to_vehss), median_ratio = median(ratio_to_vehss),
            geo_mean_ratio = exp(mean(log(ratio_to_vehss))), mean_abs_log_ratio = mean(abs(log(ratio_to_vehss))),
            rmse_pp = 100 * sqrt(mean((est - prev_vehss)^2)),
            pop_wtd_mae_pp = 100 * weighted.mean(abs(est - prev_vehss), pop_vehss),
            pearson = cor(est, prev_vehss), spearman = cor(est, prev_vehss, method = "spearman"),
            pct_in_vehss_ci = 100 * mean(in_vehss_ci), pct_ci_overlap = 100 * mean(ci_overlap),
            .groups = "drop") %>%
  mutate(method = METHOD_LABELS[as.character(method)])
message("\n--- Overall state-level agreement with VEHSS (ratio = AoU / VEHSS; pp = percentage points) ---")
print(summ %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)
write_excel_csv(summ, F_SUMMARY)

########################################################################
## Selection of combinations: most populated, highest / lowest D and E
########################################################################

eligible <- combo_summary %>% filter(n_states_DE >= MIN_STATES)
pick <- function(col, top = TRUE, n = N_SELECT) {
  x <- eligible %>% filter(is.finite(.data[[col]]))
  x <- if (top) slice_max(x, .data[[col]], n = n, with_ties = FALSE) else slice_min(x, .data[[col]], n = n, with_ties = FALSE)
  x$combo
}
sets <- list(
  pop   = pick("pop_pums", TRUE),      # most populated combinations (PUMS 40+ population)
  highD = pick("D_js", TRUE),          # largest G-vs-S divergence
  lowD  = pick("D_js", FALSE),
  highE = pick("E_js", TRUE),          # divergence spread evenly over states
  lowE  = pick("E_js", FALSE)          # divergence concentrated in a few states
)
set_labels <- c(pop = "most populated", highD = "highest D", lowD = "lowest D",
                highE = "highest E", lowE = "lowest E")
for (s in names(sets)) message(sprintf("%-15s %s", paste0(set_labels[s], ":"), paste(sets[[s]], collapse = " ; ")))

combo_summary <- combo_summary %>%
  mutate(in_pop = combo %in% sets$pop, in_highD = combo %in% sets$highD, in_lowD = combo %in% sets$lowD,
         in_highE = combo %in% sets$highE, in_lowE = combo %in% sets$lowE)

message("\n--- Combinations: population, D, E and S-RAILS agreement with VEHSS ---")
print(combo_summary %>%
        transmute(combo, n_aou, n_cases_aou, pop_pums = round(pop_pums / 1e6, 2), has_vehss, n_states_DE,
                  D = signif(D_js, 3), E = round(E_js, 3), srails_gm_ratio = round(srails_gm_ratio, 3),
                  srails_pct_in_ci = round(srails_pct_in_ci, 1)),
      n = Inf, width = Inf)

write_excel_csv(combo_state,   file.path(DATA_DIR, "vehss_cataract_combo_state_comparison.csv"))
write_excel_csv(combo_nat,     file.path(DATA_DIR, "vehss_cataract_combo_national.csv"))
write_excel_csv(combo_summary, file.path(DATA_DIR, "vehss_cataract_combo_summary.csv"))

########################################################################
## Maps (contiguous US) with REGION OUTLINES
########################################################################

us_map <- ggplot2::map_data("state") %>% rename(map_region = region)

## Region outlines: dissolve the state polygons by region with sf, then draw
## them as thick paths (works with coord_quickmap; no coord_sf needed).
region_lines <- local({
  if (!requireNamespace("sf", quietly = TRUE)) {
    message("Package 'sf' not available — region outlines are not drawn."); return(NULL)
  }
  tryCatch({
    ## planar (GEOS) geometry: the maps polygons are not valid on the sphere
    old_s2 <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(old_s2), add = TRUE)
    m  <- maps::map("state", fill = TRUE, plot = FALSE)
    sfm <- sf::st_as_sf(m)
    sfm$state  <- state.abb[match(sub(":.*$", "", sfm$ID), tolower(state.name))]
    sfm$region <- unname(state_region_map[sfm$state])
    sfm <- sf::st_make_valid(sfm[!is.na(sfm$region), ])
    ## dissolve states -> regions (aggregate() unions the geometries); a small
    ## buffer in / out closes the slivers between neighbouring state polygons
    reg <- suppressWarnings({
      r <- stats::aggregate(sfm["region"], by = list(region = sfm$region),
                            FUN = function(x) x[1], do_union = TRUE)
      sf::st_simplify(sf::st_buffer(sf::st_buffer(r, 0.02), -0.02), dTolerance = 0.005)
    })
    xy  <- sf::st_coordinates(sf::st_cast(reg, "MULTIPOLYGON"))
    tibble(long = xy[, "X"], lat = xy[, "Y"],
           group = paste(xy[, "L3"], xy[, "L2"], xy[, "L1"], sep = "_"))
  }, error = function(e) { message("Region outlines skipped: ", conditionMessage(e)); NULL })
})

poly_data <- function(tab, panel_col, value_col, panel_levels) {
  d <- tab %>%
    mutate(map_region = tolower(state.name[match(state, state.abb)]),
           panel = factor(.data[[panel_col]], levels = panel_levels),
           val   = .data[[value_col]]) %>%
    filter(!is.na(map_region), !is.na(panel))
  us_map %>%
    inner_join(d %>% select(panel, map_region, val), by = "map_region", relationship = "many-to-many")
}

choropleth <- function(tab, panel_col, value_col, panel_levels, fill_scale, title, subtitle = NULL, ncol = 2) {
  dd   <- poly_data(tab, panel_col, value_col, panel_levels)
  supp <- dplyr::filter(dd, is.na(val))
  p <- ggplot(dd, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = val), color = "grey55", linewidth = 0.08)
  if (nrow(supp) > 0) {
    if (requireNamespace("ggpattern", quietly = TRUE)) {
      p <- p + ggpattern::geom_polygon_pattern(
        data = supp, aes(long, lat, group = group),
        pattern = "stripe", pattern_angle = 45, pattern_density = 0.1,
        pattern_spacing = 0.015, pattern_size = 0.1,
        pattern_fill = "grey35", pattern_colour = NA,
        fill = "grey92", colour = "grey55", linewidth = 0.08)
    } else {
      p <- p + geom_polygon(data = supp, aes(long, lat, group = group),
                            fill = "grey75", colour = "grey55", linewidth = 0.08)
    }
  }
  if (!is.null(region_lines))          # no `panel` column -> drawn in every facet
    p <- p + geom_path(data = region_lines, aes(long, lat, group = group),
                       colour = "black", linewidth = 0.8, lineend = "round")
  p + coord_quickmap() +
    facet_wrap(~ panel, ncol = ncol) +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(strip.text = element_text(size = 9, face = "bold"),
          legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, colour = "grey35", size = 9))
}

ratio_scale <- function(vals, name) {
  L <- max(abs(vals[is.finite(vals)]), 0.1, na.rm = TRUE)
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-L, L), labels = function(x) sprintf("%.2f", 2^x),
                       na.value = "grey75", name = name)
}
prev_scale <- scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "grey75",
                                   labels = scales::percent_format(accuracy = 1), name = "Cataract\nprevalence")

slug <- function(x) gsub("[^A-Za-z0-9]+", "_", gsub("\\+", "plus", x))

## Two map figures for one combination: prevalence levels (VEHSS + 3 AoU
## estimates) and ratios (S-RAILS / VEHSS, G-RAILS / VEHSS, S-RAILS / G-RAILS)
map_combo <- function(A, S, R, tag_label = NULL) {
  st <- combo_state %>% filter(age == A, sex == S, race == R)
  cs <- combo_summary %>% filter(age == A, sex == S, race == R)
  what <- paste0(A, ", ", if (S == "Both") "both sexes" else tolower(S), ", ",
                 if (R == "All") "all races" else R)
  sub  <- paste(c(paste0("AoU n = ", format(cs$n_aou, big.mark = ","), ", Cataract cases = ",
                         format(cs$n_cases_aou, big.mark = ","),
                         if (is.finite(cs$D_js)) sprintf(";  D = %.2e, E = %.2f", cs$D_js, cs$E_js) else "",
                         ";  hatched = suppressed (<", SUPPRESS_MIN, " cases)",
                         if (!cs$has_vehss) ";  NO VEHSS EXPORT for this selection yet" else ""),
                  fallback_note), collapse = "\n")
  fn <- slug(paste(A, S, R, sep = "_"))

  lvl <- bind_rows(
    st %>% filter(method == "srails") %>% transmute(state, panel = "vehss", val = prev_vehss),
    st %>% transmute(state, panel = as.character(method), val = est)
  ) %>% mutate(panel = METHOD_LABELS[panel])
  p1 <- choropleth(lvl, "panel", "val", unname(METHOD_LABELS), prev_scale,
                   title = paste0("Cataract prevalence — ", what,
                                  if (!is.null(tag_label)) paste0("  (", tag_label, ")") else ""),
                   subtitle = sub)
  print(p1)
  ggsave(file.path(GRAPH_DIR, paste0("vehss_cataract_map_", fn, "_prevalence.png")), p1,
         width = 12, height = 8, dpi = 300, bg = "white")

  rat <- bind_rows(
    st %>% filter(method == "srails") %>% transmute(state, panel = "S-RAILS / VEHSS", val = log2_ratio),
    st %>% filter(method == "grails") %>% transmute(state, panel = "G-RAILS / VEHSS", val = log2_ratio),
    st %>% filter(method == "srails") %>% transmute(state, panel = "S-RAILS / G-RAILS", val = log2_s_g)
  )
  p2 <- choropleth(rat, "panel", "val", c("S-RAILS / VEHSS", "G-RAILS / VEHSS", "S-RAILS / G-RAILS"),
                   ratio_scale(rat$val, "ratio"), ncol = 3,
                   title = paste0("Prevalence ratios — ", what,
                                  if (!is.null(tag_label)) paste0("  (", tag_label, ")") else ""),
                   subtitle = paste0("1 = equal; blue = numerator lower, red = higher.  ", sub))
  print(p2)
  ggsave(file.path(GRAPH_DIR, paste0("vehss_cataract_map_", fn, "_ratios.png")), p2,
         width = 14, height = 4.8, dpi = 300, bg = "white")
  invisible(NULL)
}

## Overall first, then each selected set
map_combo("40+", "Both", "All", "overall")
mapped <- "40+ / Both / All"
for (s in intersect(MAP_SETS, names(sets))) {
  for (cb in setdiff(sets[[s]], mapped)) {
    parts <- strsplit(cb, " / ", fixed = TRUE)[[1]]
    map_combo(parts[1], parts[2], parts[3], set_labels[s])
    mapped <- c(mapped, cb)
  }
}
selection_tab <- bind_rows(lapply(names(sets), function(s) tibble(set = set_labels[s], combo = sets[[s]]))) %>%
  mutate(mapped = combo %in% mapped)
write_excel_csv(selection_tab, file.path(DATA_DIR, "vehss_cataract_combo_selection.csv"))

########################################################################
## Overall: scatter, state dot plot, regional comparison
########################################################################

sc <- state_long %>% filter(!suppress, !is.na(prev_vehss)) %>%
  mutate(panel = factor(METHOD_LABELS[as.character(method)], levels = METHOD_LABELS[names(METHODS)]))
p_sc <- ggplot(sc, aes(prev_vehss, est)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, colour = "grey70", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lb_vehss, xmax = ub_vehss), width = 0, colour = "grey70", linewidth = 0.3, na.rm = TRUE) +
  geom_point(aes(size = pop_vehss, colour = method), alpha = 0.85) +
  scale_colour_manual(values = METHOD_COLS[names(METHODS)], guide = "none") +
  scale_size_continuous(name = "Medicare 40+\nbeneficiaries (VEHSS)", range = c(1, 6),
                        labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "VEHSS modeled crude prevalence (2019)", y = "AoU estimate",
       title = paste0("State Cataract prevalence, adults ", AGE_MIN, "+: AoU vs VEHSS"), subtitle = fallback_note) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_sc <- p_sc + ggrepel::geom_text_repel(aes(label = state), size = 2.6, colour = "grey20",
                                          max.overlaps = 20, seed = 1, segment.size = 0.2)
}
print(p_sc)
ggsave(file.path(GRAPH_DIR, "vehss_cataract_scatter_state.png"), p_sc, width = 13, height = 5, dpi = 300, bg = "white")

dp_order <- vehss_state %>% filter(!is.na(prev_vehss)) %>% arrange(prev_vehss) %>% pull(state)
dp_off   <- c(unweighted = 0.22, grails = 0, srails = -0.22)
dp <- state_long %>% filter(!suppress, state %in% dp_order) %>%
  mutate(y = match(state, dp_order) + dp_off[as.character(method)],
         method_lab = factor(METHOD_LABELS[as.character(method)], levels = METHOD_LABELS))
dp_v <- vehss_state %>% filter(state %in% dp_order) %>% mutate(y = match(state, dp_order),
                               method_lab = factor(METHOD_LABELS["vehss"], levels = METHOD_LABELS))
p_dp <- ggplot() +
  geom_linerange(data = dp_v, aes(y = y, xmin = lb_vehss, xmax = ub_vehss), colour = "grey80", linewidth = 2.2, na.rm = TRUE) +
  geom_point(data = dp_v, aes(y = y, x = prev_vehss, colour = method_lab), shape = 18, size = 3) +
  geom_point(data = dp, aes(y = y, x = est, colour = method_lab), size = 1.6) +
  scale_colour_manual(values = setNames(METHOD_COLS, METHOD_LABELS), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(breaks = seq_along(dp_order), labels = dp_order, expand = expansion(add = 0.6)) +
  labs(x = paste0("Cataract prevalence, adults ", AGE_MIN, "+"), y = NULL,
       title = "State Cataract prevalence: AoU estimates against the VEHSS 95% interval",
       subtitle = paste(c(paste0("Grey band = VEHSS 95% CI; suppressed states (<", SUPPRESS_MIN,
                                 " AoU cases) show VEHSS only"), fallback_note), collapse = "\n")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
print(p_dp)
ggsave(file.path(GRAPH_DIR, "vehss_cataract_dotplot_state.png"), p_dp, width = 8, height = 11, dpi = 300, bg = "white")

reg_plot <- bind_rows(
  reg %>% transmute(area = as.character(region), method = as.character(method), est, lb, ub),
  reg %>% distinct(region, prev_vehss, lb_vehss, ub_vehss) %>%
    transmute(area = as.character(region), method = "vehss", est = prev_vehss, lb = lb_vehss, ub = ub_vehss),
  nat %>% transmute(area = "National", method = as.character(method), est, lb, ub),
  vehss_us %>% transmute(area = "National", method = "vehss", est = prev_vehss, lb = lb_vehss, ub = ub_vehss)
) %>%
  mutate(area = factor(area, levels = c(REGION_LEVELS, "National")),
         method_lab = factor(METHOD_LABELS[method], levels = METHOD_LABELS),
         lb = coalesce(lb, est), ub = coalesce(ub, est))   # regional VEHSS has no CI
p_reg <- ggplot(reg_plot, aes(x = area, y = est, colour = method_lab)) +
  geom_pointrange(aes(ymin = lb, ymax = ub), position = position_dodge(width = 0.6), size = 0.4) +
  scale_colour_manual(values = setNames(METHOD_COLS, METHOD_LABELS), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = paste0("Cataract prevalence, adults ", AGE_MIN, "+"),
       title = "Regional and national Cataract prevalence: AoU vs VEHSS",
       subtitle = paste(c("Regional VEHSS = state estimates aggregated by implied 40+ population (no CI); AoU bars = 95% CI",
                          fallback_note), collapse = "\n")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
print(p_reg)
ggsave(file.path(GRAPH_DIR, "vehss_cataract_region_comparison.png"), p_reg, width = 9, height = 5.5, dpi = 300, bg = "white")

########################################################################
## Combinations: S-RAILS vs VEHSS scatter for the mapped combinations, and
## the D-vs-E scatter of all combinations
########################################################################

sc_c <- combo_state %>%
  mutate(combo = paste(age, sex, race, sep = " / ")) %>%
  filter(combo %in% mapped, method == "srails", !suppress, !is.na(prev_vehss)) %>%
  mutate(combo = factor(combo, levels = mapped))
if (nrow(sc_c)) {
  p_scc <- ggplot(sc_c, aes(prev_vehss, est)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lb_vehss, xmax = ub_vehss), width = 0, colour = "grey75", linewidth = 0.3, na.rm = TRUE) +
    geom_point(aes(size = pop_vehss), colour = METHOD_COLS["srails"], alpha = 0.85) +
    scale_size_continuous(name = "population\n(VEHSS)", range = c(1, 5),
                          labels = scales::label_number(scale = 1e-6, suffix = "M")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    facet_wrap(~ combo, scales = "free", labeller = label_wrap_gen(28)) +
    labs(x = "VEHSS modeled crude prevalence (2019)", y = "AoU S-RAILS (region) estimate",
         title = "State Cataract prevalence by Age x Sex x Race: S-RAILS vs VEHSS", subtitle = fallback_note) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  if (requireNamespace("ggrepel", quietly = TRUE))
    p_scc <- p_scc + ggrepel::geom_text_repel(aes(label = state), size = 2.2, colour = "grey20",
                                              max.overlaps = 12, seed = 1, segment.size = 0.2)
  print(p_scc)
  ggsave(file.path(GRAPH_DIR, "vehss_cataract_scatter_combos_srails_vs_vehss.png"), p_scc,
         width = 13, height = 9, dpi = 300, bg = "white")
}

## D vs E: every combination with a finite D and E; highlight the selection
## sets by colour, the most populated by a ring, the overall as a diamond.
SET_COLOURS <- c("highest D" = "#C0392B", "lowest D" = "#2C7FB8", "highest E" = "#1B7837",
                 "lowest E" = "#E08214", "other" = "grey82")
de <- combo_summary %>%
  filter(is.finite(D_js), is.finite(E_js), n_states_DE >= MIN_STATES) %>%
  mutate(set = case_when(in_highD ~ "highest D", in_lowE ~ "lowest E", in_highE ~ "highest E",
                         in_lowD ~ "lowest D", TRUE ~ "other"),
         set = factor(set, levels = names(SET_COLOURS)),
         populated = ifelse(in_pop, "most populated", "other"),
         shape = ifelse(is_overall, "overall (40+, both sexes, all races)", "combination"))
n_hidden <- nrow(combo_summary) - nrow(de)     # too few states, or G-RAILS == S-RAILS everywhere
lab_de <- de %>% filter(set != "other" | in_pop | is_overall)
p_de <- ggplot(de, aes(D_js, E_js)) +
  geom_point(aes(fill = set, size = pop_pums, colour = populated, shape = shape), alpha = 0.85, stroke = 0.7) +
  scale_fill_manual(values = SET_COLOURS, name = "selection set", breaks = setdiff(names(SET_COLOURS), "other")) +
  scale_colour_manual(values = c("most populated" = "grey10", "other" = "#00000000"),   # transparent ring, not NA (NA drops the point)
                      breaks = "most populated", name = NULL) +
  scale_shape_manual(values = c("combination" = 21, "overall (40+, both sexes, all races)" = 23), name = NULL) +
  scale_size_continuous(name = "40+ population\n(PUMS)", range = c(1.5, 8),
                        labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_log10(labels = function(v) format(v, scientific = TRUE, digits = 2)) +
  labs(x = expression(D[j] ~ "(population-weighted G-RAILS vs S-RAILS divergence, log scale)"),
       y = expression(E[j] ~ "(dispersion across states)"),
       title = "Divergence D and dispersion E of the G-RAILS / S-RAILS difference, by Age x Sex x Race",
       subtitle = paste(c("Right = the two weightings disagree more; low = the disagreement sits in a few states.",
                          paste0(nrow(de), " combinations shown; ", n_hidden, " with fewer than ", MIN_STATES,
                                 " un-suppressed states or D = 0 are not shown."),
                          fallback_note), collapse = "\n")) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 4, colour = NA, shape = 21)),
         colour = guide_legend(order = 2, override.aes = list(size = 4, fill = "grey82", stroke = 1, shape = 21)),
         shape = guide_legend(order = 3, override.aes = list(size = 4, fill = "grey82")),
         size = guide_legend(order = 4, override.aes = list(fill = "grey70", colour = NA, shape = 21))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey35"), legend.key = element_blank())
if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab_de))
  p_de <- p_de + ggrepel::geom_text_repel(data = lab_de, aes(label = combo), size = 2.8, colour = "grey15",
                                          box.padding = 0.5, point.padding = 0.3, min.segment.length = 0,
                                          segment.size = 0.3, segment.colour = "grey55", max.overlaps = Inf, seed = 1)
print(p_de)
ggsave(file.path(GRAPH_DIR, "vehss_cataract_DE_scatter_combinations.png"), p_de, width = 11, height = 7.5, dpi = 300, bg = "white")

message("\nDone. Tables -> ", DATA_DIR, "/vehss_*.csv ; figures -> ", GRAPH_DIR, "/vehss_*.png",
        if (length(fallback_regions)) paste0("\nREMINDER: ", fallback_note) else "")
