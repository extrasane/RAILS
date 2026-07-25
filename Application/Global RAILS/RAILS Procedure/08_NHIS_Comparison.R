#### 08_NHIS_Comparison.R — AoU vs NHIS 2020 health-outcome comparison
#### Run inside the AoU Researcher Workbench.
####
#### Produces the manuscript table comparing three dynamic health outcomes
#### between NHIS 2020 (design-weighted) and All of Us under three schemes:
####   Unweighted        — raw AoU proportions
####   Recalibrate       — RAILS weights (d_rails)
####   Double-Weighting  — RAILS weights divided by an item-response
####                       propensity, correcting for question non-response
####
#### NOTE (updated Workbench): the person/survey/state queries no longer
#### restrict to has_ehr_data = 1, and the survey pull adds concept IDs
#### 1585386 (insurance) and 1585711 (overall health).
####
#### Workbench 2.0: queries download straight into R via bq_project_query +
#### bq_table_download. The old bq_table_save route is avoided because it
#### must stage results in a GCS bucket, and WORKSPACE_BUCKET can point at a
#### bucket that does not exist (e.g. "gs://cloned-mybucket-..." inherited by
#### a cloned workspace), which fails with "Not found: URI ... [notFound]".

library(tidyverse)
library(bigrquery)
library(lubridate)
library(survey)
library(xtable)

########################################################################
## Workbench 2.0 / BigQuery setup
########################################################################

## Set ONE of these before running (Resources panel gives the resource id):
##   Sys.setenv(AOU_CDR_RESOURCE_ID = "my-controlled-tier-cdr")
##   Sys.setenv(AOU_CDR_DATASET     = "project.dataset")
AOU_CDR_RESOURCE_ID <- Sys.getenv("AOU_CDR_RESOURCE_ID", unset = "")
AOU_CDR_DATASET     <- Sys.getenv("AOU_CDR_DATASET",     unset = "")
BQ_BILLING_PROJECT  <- Sys.getenv("BQ_BILLING_PROJECT",  unset = "")

is_nonempty <- function(x) !is.na(x) && nzchar(trimws(x))

run_cmd <- function(cmd, args) {
  out <- tryCatch(system2(cmd, args = args, stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
  out[nzchar(out)]
}

resolve_wb_resource <- function(resource_id) {
  if (!is_nonempty(resource_id)) return(NA_character_)
  for (args in list(c("resolve", paste0("--id=", resource_id)),
                    c("resource", "resolve", paste0("--id=", resource_id)))) {
    out <- run_cmd("wb", args)
    hit <- out[grepl("\\.", out)]
    if (length(hit) > 0) return(hit[1])
  }
  stop("Could not resolve Workbench resource id: ", resource_id)
}

get_cdr_dataset <- function() {
  if (is_nonempty(AOU_CDR_DATASET)) return(AOU_CDR_DATASET)
  legacy <- Sys.getenv("WORKSPACE_CDR", unset = "")
  if (is_nonempty(legacy)) return(legacy)
  if (is_nonempty(AOU_CDR_RESOURCE_ID)) return(resolve_wb_resource(AOU_CDR_RESOURCE_ID))
  stop("No AoU CDR dataset found. Set AOU_CDR_RESOURCE_ID or AOU_CDR_DATASET.")
}

parse_bq_dataset <- function(dataset_path) {
  parts <- strsplit(gsub("`", "", dataset_path), "\\.")[[1]]
  if (length(parts) != 2)
    stop("Expected 'project.dataset', got: ", dataset_path)
  list(project = parts[1], dataset = parts[2], path = paste(parts, collapse = "."))
}

get_billing_project <- function(default_project) {
  cand <- c(BQ_BILLING_PROJECT,
            Sys.getenv("GOOGLE_PROJECT",       unset = ""),
            Sys.getenv("GOOGLE_CLOUD_PROJECT", unset = ""),
            Sys.getenv("GCLOUD_PROJECT",       unset = ""))
  cand <- cand[nzchar(cand)]
  if (length(cand) > 0) return(cand[1])
  default_project
}

cdr             <- parse_bq_dataset(get_cdr_dataset())
billing_project <- get_billing_project(cdr$project)
message("Using AoU CDR dataset: ", cdr$path)
message("Using BigQuery billing project: ", billing_project)

fq_table <- function(table_name) paste0("`", cdr$path, ".", table_name, "`")

## bigint = "numeric" returns person_id as a plain double, matching what
## read_csv() produces for the weights file. With the "integer64" default the
## join fails: "Can't join `x$person_id` with `y$person_id` due to
## incompatible types". AoU person_ids are ~9 digits, far below 2^53, so
## double representation is exact.
run_bq <- function(sql, bigint = "numeric") {
  message("Submitting BigQuery job...")
  job <- bq_project_query(x = billing_project, query = sql, use_legacy_sql = FALSE)
  bq_table_download(job, bigint = bigint, quiet = TRUE)
}

########################################################################
## BigQuery extraction
########################################################################

cohort_sql <- paste0("SELECT DISTINCT person_id FROM ", fq_table("cb_search_person"))

dataset_person_sql <- paste0("
  SELECT
    person.person_id,
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
  WHERE person.person_id IN (", cohort_sql, ")")

## Survey concept IDs: income, homeown, edu, careplace, insurance, health
dataset_survey_sql <- paste0("
  SELECT answer.person_id, answer.question_concept_id, answer.question, answer.answer
  FROM ", fq_table("ds_survey"), " answer
  WHERE question_concept_id IN (1585370, 1585375, 1585899, 1585940,
                                1585386, 1585711, 43530593)
    AND answer.person_id IN (", cohort_sql, ")")

dataset_state_sql <- paste0("
  SELECT person.person_id, person.state_of_residence_source_value AS state
  FROM ", fq_table("person_ext"), " person
  WHERE person.person_id IN (", cohort_sql, ")")

## person_id is forced to double everywhere. BigQuery may return it as
## <integer64>, which will not join against the <double> that read_csv gives
## for the weights file. Use as.numeric() (dispatches to bit64's method) —
## never unclass(), which reinterprets the bit pattern and yields garbage.
person_df <- run_bq(dataset_person_sql) %>%
  mutate(person_id = as.numeric(person_id),
         across(c(gender, race, ethnicity, sex_at_birth), as.character))
survey_df <- run_bq(dataset_survey_sql) %>%
  mutate(person_id = as.numeric(person_id),
         across(c(question, answer), as.character))
state_df  <- run_bq(dataset_state_sql) %>%
  mutate(person_id = as.numeric(person_id),
         state = as.character(state))

message("Rows downloaded: person=", nrow(person_df),
        ", survey=", nrow(survey_df), ", state=", nrow(state_df))

########################################################################
## Clean and recode
########################################################################

survey_wide <- survey_df %>%
  select(!question_concept_id) %>%
  filter(question != "The Basics: Sexual Orientation") %>%
  pivot_wider(names_from = question, values_from = answer, values_fn = first)

dt0 <- person_df %>%
  full_join(survey_wide, by = "person_id") %>%
  full_join(state_df,    by = "person_id") %>%
  rename(income        = `Income: Annual Income`,
         homeown       = `Home Own: Current Home Own`,
         edu           = `Education Level: Highest Grade`,
         careplace     = `Health Advice: What Kind Of Place`,
         insur         = `Insurance: Health Insurance`,
         health_status = `Overall Health: General Health`)

dt0 <- dt0 %>%
  mutate(date_of_birth = as.Date(date_of_birth),
         age = round(interval(date_of_birth, as.Date("2024-08-01")) /
                       duration(num = 1, units = "years"), 0)) %>%
  select(!date_of_birth) %>%
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
    ethnicity %in% c("PMI: Prefer Not To Answer", "PMI: Skip")            ~ NA_character_)) %>%
  mutate(sex = case_when(
    sex_at_birth %in% c("I prefer not to answer", "Intersex", "No matching concept",
                        "None", "PMI: Skip") ~ NA_character_,
    TRUE                                     ~ sex_at_birth)) %>%
  select(!sex_at_birth) %>%
  mutate(income = case_when(
    income %in% c("Annual Income: 100k 150k", "Annual Income: 150k 200k",
                  "Annual Income: more 200k")                  ~ ">100k",
    income %in% c("Annual Income: 75k 100k")                   ~ "75k-100k",
    income %in% c("Annual Income: 50k 75k")                    ~ "50k-75k",
    income %in% c("Annual Income: 35k 50k")                    ~ "35k-50k",
    income %in% c("Annual Income: 25k 35k", "Annual Income: 10k 25k",
                  "Annual Income: less 10k")                   ~ "<35k",
    TRUE                                                       ~ NA_character_)) %>%
  mutate(homeown = case_when(
    homeown %in% c("Current Home Own: Other Arrangement", "PMI: Dont Know") ~ "Other",
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
  mutate(careplace = case_when(
    careplace == "What Kind Of Place: Doctors Office"          ~ "Doctors office",
    careplace == "What Kind Of Place: Emergency Room"          ~ "Emergency",
    careplace == "What Kind Of Place: No One Place Most Often" ~ "None",
    careplace == "What Kind Of Place: Urgent Care"             ~ "Urgent",
    careplace == "What Kind Of Place: Some Other Place"        ~ "Others",
    TRUE                                                       ~ NA_character_)) %>%
  mutate(insur = case_when(
    insur == "Health Insurance: No"  ~ "Not covered",
    insur == "Health Insurance: Yes" ~ "Covered",
    TRUE                             ~ NA_character_)) %>%
  ## "Overall Health: Excellent" -> "Excellent"
  mutate(health_status = str_replace_all(health_status, "^[^:]*: ", ""),
         health_status = ifelse(health_status == "Skip", NA_character_, health_status))

dt00 <- dt0 %>%
  mutate(state = ifelse(str_detect(state, "^PII State: "),
                        str_extract(state, "(?<=PII State: )\\w{2}"), NA))

########################################################################
## Attach RAILS weights and apply factor levels
########################################################################

## RAILS weights — read from the local disk.
##
## ## BUCKET ROUTE (annotated out): fetches the file from GCS if it is not
## ## already on the Jupyter disk. WORKSPACE_BUCKET can be stale in a cloned
## ## workspace and point at a bucket that does not exist, so set
## ## WEIGHTS_BUCKET or OUTPUT_BUCKET_RESOURCE_ID before enabling this.
## WEIGHTS_FILE <- "all_dt_rails.csv"
##
## if (!file.exists(WEIGHTS_FILE)) {
##   bucket <- Sys.getenv("WEIGHTS_BUCKET", unset = "")
##   if (!is_nonempty(bucket) && is_nonempty(Sys.getenv("OUTPUT_BUCKET_RESOURCE_ID", unset = "")))
##     bucket <- resolve_wb_resource(Sys.getenv("OUTPUT_BUCKET_RESOURCE_ID"))
##   if (!is_nonempty(bucket)) bucket <- Sys.getenv("WORKSPACE_BUCKET", unset = "")
##   if (!is_nonempty(bucket))
##     stop("Could not locate ", WEIGHTS_FILE, ". Upload it to the working directory, ",
##          "or set WEIGHTS_BUCKET / OUTPUT_BUCKET_RESOURCE_ID.")
##   message("Fetching ", WEIGHTS_FILE, " from ", bucket)
##   run_cmd("gsutil", c("cp", paste0(sub("/$", "", bucket), "/data/", WEIGHTS_FILE), "."))
## }
## raking.wts <- read_csv(WEIGHTS_FILE, show_col_types = FALSE)

raking.wts <- read_csv("../data/global_rails_weights.csv", show_col_types = FALSE)

## Population total = the sum of the RAILS weights themselves, so the scale
## used downstream matches the weights actually being applied.
nsiz <- sum(raking.wts$w_rails)

## Keep only the key and the RAILS weight; the other benchmark columns
## (w_unweighted, w_cal1, ...) are not used in this comparison.
raking.wts <- raking.wts %>%
  select(person_id, w_rails)

south     <- c("AL","AR","FL","GA","KY","LA","MS","NC","SC","TN","TX","VA","WV")
midwest   <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI")
northeast <- c("CT","DE","ME","MD","MA","NH","NJ","NY","PA","RI","VT")
west      <- c("AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY")
state_region_map <- c(
  setNames(rep("South",     length(south)),     south),
  setNames(rep("Midwest",   length(midwest)),   midwest),
  setNames(rep("Northeast", length(northeast)), northeast),
  setNames(rep("West",      length(west)),      west))

dt000 <- dt00 %>%
  mutate(race_eth = case_when(
    ethnicity == "Yes"                       ~ "Hispanic",
    race      == "Asian"                     ~ "NH Asian",
    race      == "Black or African American" ~ "NH Black",
    race      == "White"                     ~ "NH White",
    TRUE                                     ~ "Others")) %>%
  mutate(agegroup = case_when(
    age <= 24            ~ "18-24",
    age > 24 & age <= 44 ~ "25-44",
    age > 44 & age <= 64 ~ "45-64",
    age > 64 & age <= 74 ~ "65-74",
    TRUE                 ~ "75+")) %>%
  mutate(sex           = factor(sex,      levels = c("Female","Male")),
         race_eth      = factor(race_eth, levels = c("Hispanic","NH Asian","NH Black","NH White","Others")),
         income        = factor(income,   levels = c("<35k","35k-50k","50k-75k","75k-100k",">100k")),
         agegroup      = factor(agegroup, levels = c("18-24","25-44","45-64","65-74","75+")),
         edu           = factor(edu,      levels = c("Less than highschool","Some highschool",
                                                     "Highschool graduate","Some college",
                                                     "College graduate or advanced")),
         homeown       = factor(ifelse(homeown == "Other", "Others", homeown),
                                levels = c("Own","Rent","Others")),
         health_status = factor(health_status, levels = c("Excellent","Very Good","Good","Fair","Poor")),
         careplace     = factor(careplace,     levels = c("Doctors office","Urgent","Emergency","Others","None")),
         insur         = factor(insur,         levels = c("Not covered","Covered")),
         region        = factor(state_region_map[state],
                                levels = c("Northeast","Midwest","South","West"))) %>%
  ## person_id is already double on both sides (coerced at download and by
  ## read_csv), so this joins cleanly.
  full_join(raking.wts, by = "person_id") %>%
  rename(weight = w_rails) %>%
  select(-c(gender, ethnicity)) %>%
  mutate(COM = ifelse(is.na(weight), 0, 1))

## nsiz is set above from sum(raking.wts$w_rails)
message("nsiz (sum of RAILS weights): ", format(nsiz, big.mark = ","))
message("Participants with a RAILS weight: ", sum(dt000$COM), " of ", nrow(dt000))

########################################################################
## NHIS 2020 reference
########################################################################

zip_url   <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NHIS/2020/adult20csv.zip"
temp      <- tempfile(); download.file(zip_url, temp)
unzip_dir <- tempdir(); unzip(temp, exdir = unzip_dir)
df        <- read.csv(list.files(unzip_dir, pattern = "20.csv$", full.names = TRUE))

dt_nhis <- df %>%
  select(NOTCOV_A, USPLKIND_A, PHSTAT_A, WTFA_A) %>%
  rename(insur = NOTCOV_A, careplace = USPLKIND_A,
         health_status = PHSTAT_A, weight = WTFA_A) %>%
  mutate(insur = case_when(insur == 1 ~ "Not covered",
                           insur == 2 ~ "Covered",
                           TRUE       ~ NA_character_)) %>%
  mutate(careplace = case_when(careplace == 1     ~ "Doctors office",
                               careplace == 2     ~ "Urgent",
                               careplace == 3     ~ "Emergency",
                               careplace %in% 4:5 ~ "Others",
                               careplace == 6     ~ "None",
                               TRUE               ~ NA_character_)) %>%
  mutate(health_status = case_when(health_status == 1 ~ "Excellent",
                                   health_status == 2 ~ "Very Good",
                                   health_status == 3 ~ "Good",
                                   health_status == 4 ~ "Fair",
                                   health_status == 5 ~ "Poor",
                                   TRUE               ~ NA_character_)) %>%
  mutate(health_status = factor(health_status, levels = c("Excellent","Very Good","Good","Fair","Poor")),
         careplace     = factor(careplace,     levels = c("Doctors office","Urgent","Emergency","Others","None")),
         insur         = factor(insur,         levels = c("Not covered","Covered")))

########################################################################
## Estimation helpers
########################################################################

logit <- function(x) log(x / (1 - x))
expit <- function(x) 1 / (1 + exp(-x))

## Outcomes shown in the table (in table order)
names_univar <- c("careplace", "health_status", "insur")

## CI level. NOTE: qnorm(0.95) gives a 90% interval; the manuscript table is
## labelled "95% CI", which requires qnorm(0.975).
Z <- qnorm(0.975)

## Logit-scale CI for a proportion, as defined in the paper: the size used in
## the variance is the estimate's own sum of weights (the unweighted count
## when all weights are 1). Vectorized — p and size are vectors with one
## element per level, and each function returns a vector of the same length.
ci_se <- function(p, size) sqrt(1 / (size * p * (1 - p)))
ci_lb <- function(p, size) expit(logit(p) - Z * ci_se(p, size))
ci_ub <- function(p, size) expit(logit(p) + Z * ci_se(p, size))

## Per-level summary for one variable under one weight column
summarize_var <- function(df, var, wt, label, scale = 1e7) {
  d <- df %>% filter(!is.na(.data[[var]]))
  n_tot <- nrow(d)
  w_tot <- sum(d[[wt]], na.rm = TRUE)
  d %>%
    group_by(Level = .data[[var]]) %>%
    summarise(n = n(), wsum = sum(.data[[wt]], na.rm = TRUE), .groups = "drop") %>%
    mutate(Variable = var,
           prev = wsum / w_tot,
           ## Variance size = the level's own sum of weights (equals the
           ## unweighted count when the weight column is all 1s)
           lb   = ci_lb(prev, wsum),
           ub   = ci_ub(prev, wsum),
           size = wsum / scale) %>%
    select(Variable, Level, size, prev, lb, ub) %>%
    rename_with(~paste0(label, "_", .x), c(size, prev, lb, ub))
}

########################################################################
## (1) NHIS — design-weighted
########################################################################

margin_nhis <- map_dfr(names_univar,
                       ~summarize_var(dt_nhis, .x, "weight", "NHIS"))

########################################################################
## (2) AoU unweighted, and (3) AoU RAILS-weighted ("Recalibrate")
########################################################################

## Weights are rescaled to nsiz PER OUTCOME, after filtering to that
## outcome's respondents — matching the original workflow. (Rescaling once
## globally leaves the prevalences unchanged, since the scale cancels in
## wsum / w_tot, but changes the reported Size column.)
aou_unw <- map_dfr(names_univar, function(var) {
  dt000 %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(wt_one = 1) %>%
    summarize_var(var, "wt_one", "UNW", scale = 1e3)
})

aou_w <- map_dfr(names_univar, function(var) {
  dt000 %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(weight_scaled = weight * nsiz / sum(weight, na.rm = TRUE)) %>%
    summarize_var(var, "weight_scaled", "W")
})

########################################################################
## (4) AoU double-weighted ("Double-Weighting")
##   RAILS weight / P(answered the item), the response propensity modelled
##   on the demographic covariates. Corrects for item non-response on top
##   of the participation correction RAILS already applies.
########################################################################

ps_covars <- c("age", "sex", "race_eth", "region", "homeown", "edu", "income")

aou_rw <- map_dfr(names_univar, function(var) {
  temp <- dt000 %>%
    select(person_id, all_of(var), all_of(ps_covars), COM, weight) %>%
    mutate(index = ifelse(is.na(.data[[var]]), 0, 1))

  f     <- as.formula(paste("index ~", paste(ps_covars, collapse = "+")))
  model <- glm(f, data = temp, family = binomial)

  temp %>%
    filter(if_all(all_of(ps_covars), ~!is.na(.x))) %>%
    mutate(ps = predict(model, newdata = ., type = "response"),
           weight = weight / ps,
           weight = weight * nsiz / sum(weight, na.rm = TRUE)) %>%
    filter(!is.na(.data[[var]])) %>%
    summarize_var(var, "weight", "RW")
})

########################################################################
## Assemble the comparison table
########################################################################

out <- margin_nhis %>%
  left_join(aou_unw, by = c("Variable", "Level")) %>%
  left_join(aou_w,   by = c("Variable", "Level")) %>%
  left_join(aou_rw,  by = c("Variable", "Level")) %>%
  mutate(
    ## Diff = NHIS - AoU, in percentage points
    diff_unw = 100 * (NHIS_prev - UNW_prev),
    diff_w   = 100 * (NHIS_prev - W_prev),
    diff_rw  = 100 * (NHIS_prev - RW_prev)
  )

## Percentages for display
tab <- out %>%
  mutate(across(matches("_prev$|_lb$|_ub$"), ~ .x * 100)) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  select(Variable, Level,
         NHIS_size, NHIS_prev, NHIS_lb, NHIS_ub,
         UNW_size,  UNW_prev,  UNW_lb,  UNW_ub,  diff_unw,
         W_size,    W_prev,    W_lb,    W_ub,    diff_w,
         RW_size,   RW_prev,   RW_lb,   RW_ub,   diff_rw)

print(tab, n = Inf)

## Total absolute difference — the "Total Absolute Diff" row of each block
totals <- tab %>%
  group_by(Variable) %>%
  summarise(across(starts_with("diff"), ~sum(abs(.x))), .groups = "drop")
print(totals)

## Overall, across all outcomes
print(tab %>% summarise(across(starts_with("diff"), ~sum(abs(.x)))))

########################################################################
## LaTeX output
##   Percent columns are pre-formatted as "xx.xx (lb - ub)" so they drop
##   straight into the manuscript table's Percent (95% CI) cells.
########################################################################

fmt_ci <- function(p, lb, ub) sprintf("%.2f \\, (%.2f - %.2f)", p, lb, ub)

tex <- tab %>%
  transmute(
    Level,
    NHIS_size, NHIS_pct = fmt_ci(NHIS_prev, NHIS_lb, NHIS_ub),
    UNW_size,  UNW_pct  = fmt_ci(UNW_prev,  UNW_lb,  UNW_ub),  diff_unw,
    W_size,    W_pct    = fmt_ci(W_prev,    W_lb,    W_ub),    diff_w,
    RW_size,   RW_pct   = fmt_ci(RW_prev,   RW_lb,   RW_ub),   diff_rw)

print(xtable(tex), type = "latex", include.rownames = FALSE,
      sanitize.text.function = identity)

write_excel_csv(tab, "nhis_aou_comparison.csv")
# Optional upload (needs a valid bucket — see the WEIGHTS_BUCKET note above):
# run_cmd("gsutil", c("cp", "./nhis_aou_comparison.csv",
#                     paste0(sub("/$", "", Sys.getenv("WEIGHTS_BUCKET")), "/data/")))
