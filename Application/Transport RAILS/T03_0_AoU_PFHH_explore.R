#### T03_0_AoU_PFHH_explore.R — Part 3, step 0: how are the health-history
#### survey items stored in this CDR?
#### Run inside the AoU Researcher Workbench (Controlled Tier CDR attached).
####
#### Purpose: fix the "AoU survey" column of the Part 3 outcome crosswalk
#### (ROADMAP_Part3.md, section 4). Nothing is analysed here; the script only
#### LISTS the survey(s), questions and answers of the Personal and Family
#### Health History (and the older Personal / Family Medical History) with the
#### number of distinct participants, and flags the items that match the
#### crosswalk concepts (hypertension, cholesterol, diabetes, CHD, MI, angina,
#### stroke, asthma, COPD, cancer, arthritis, depression, anxiety, kidney
#### disease, hepatitis).
####
#### Outputs (../data/):
####   T03_0_ds_survey_columns.csv        column names of ds_survey (structure check)
####   T03_0_survey_names.csv             every survey name with participants and date range
####   T03_0_health_history_items.csv     every question x answer of the health-history survey(s)
####   T03_0_crosswalk_candidates.csv     the subset matching the crosswalk concepts
#### Counts below SUPPRESS_MIN are replaced by NA (AoU disclosure policy), so the
#### csvs may be shared; please send T03_0_crosswalk_candidates.csv (and, if
#### convenient, T03_0_survey_names.csv) back for the crosswalk.

suppressPackageStartupMessages({
  library(tidyverse)
  library(bigrquery)
  library(jsonlite)
})

DATA_DIR     <- "../data"
SUPPRESS_MIN <- 20
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

## Surveys considered health history (case-insensitive regex on ds_survey.survey)
HISTORY_SURVEY_REGEX <- "health history|medical history"

## Crosswalk concepts -> keyword regex on question OR answer text (lower case)
CONCEPTS <- c(
  hypertension     = "hypertension|high blood pressure",
  high_cholesterol = "cholesterol|hyperlipid",
  diabetes         = "diabetes",
  chd              = "coronary",
  mi               = "heart attack|myocardial",
  angina           = "angina",
  stroke           = "stroke",
  asthma           = "asthma",
  copd             = "copd|emphysema|chronic bronchitis|chronic obstructive",
  cancer           = "cancer",
  arthritis        = "arthritis",
  depression       = "depress",
  anxiety          = "anxiety",
  kidney           = "kidney",
  hepatitis        = "hepatitis"
)

########################################################################
## Workbench 2.0 CDR resolution (identical to the VEHSS_comparison*.R scripts)
########################################################################

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
wb_json <- function(args) {
  out <- run_cmd("wb", args)
  if (!length(out)) return(NULL)
  txt   <- paste(out, collapse = "\n")
  start <- regexpr("[\\[{]", txt)
  if (start < 0) return(NULL)
  tryCatch(fromJSON(substr(txt, start, nchar(txt)), simplifyVector = FALSE), error = function(e) NULL)
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
list_bq_resources <- function() {
  res <- wb_json(c("resource", "list", "--type=BQ_DATASET", "--format=JSON"))
  if (is.null(res)) return(NULL)
  if (!is.null(names(res)) && !is.null(res$id)) res <- list(res)
  tab <- bind_rows(lapply(res, function(r) tibble(
    id          = as.character(r$id %||% r$name %||% NA_character_),
    stewardship = as.character(r$stewardshipType %||% r$stewardship %||% NA_character_),
    description = as.character(r$description %||% ""),
    projectId   = as.character(r$projectId %||% r$resourceAttributes$gcpBqDataset$projectId %||% NA_character_),
    datasetId   = as.character(r$datasetId %||% r$resourceAttributes$gcpBqDataset$datasetId %||% NA_character_))))
  if (nrow(tab) == 0) NULL else tab
}
get_cdr_dataset <- function() {
  if (is_nonempty(AOU_CDR_DATASET))     return(AOU_CDR_DATASET)
  if (is_nonempty(AOU_CDR_RESOURCE_ID)) return(resolve_wb_resource(AOU_CDR_RESOURCE_ID))
  tab <- list_bq_resources()
  if (!is.null(tab)) {
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
      if (!is.na(cand$projectId) && !is.na(cand$datasetId)) return(paste0(cand$projectId, ".", cand$datasetId))
      return(resolve_wb_resource(cand$id))
    }
    stop("Several BigQuery dataset resources could be the CDR: ", paste(cand$id, collapse = ", "),
         ".\nPick one with  Sys.setenv(AOU_CDR_RESOURCE_ID = '<resource id>')  and re-run.")
  }
  legacy <- Sys.getenv("WORKSPACE_CDR", unset = "")
  if (is_nonempty(legacy)) return(legacy)
  stop("No AoU CDR dataset found. Set Sys.setenv(AOU_CDR_DATASET = 'project.dataset').")
}
parse_bq_dataset <- function(dataset_path) {
  parts <- strsplit(gsub("`", "", trimws(dataset_path)), "\\.")[[1]]
  if (length(parts) != 2) stop("Expected 'project.dataset', got: ", dataset_path)
  list(project = parts[1], dataset = parts[2], path = paste(parts, collapse = "."))
}
get_billing_project <- function(default_project) {
  cand <- c(BQ_BILLING_PROJECT, Sys.getenv("GOOGLE_CLOUD_PROJECT", unset = ""),
            Sys.getenv("GOOGLE_PROJECT", unset = ""), Sys.getenv("GCLOUD_PROJECT", unset = ""))
  cand <- cand[nzchar(cand)]
  if (length(cand) > 0) return(cand[1])
  ws <- wb_json(c("workspace", "describe", "--format=JSON"))
  wp <- ws$googleProjectId %||% ws$gcpProjectId %||% ws$projectId %||% NULL
  if (is_nonempty(wp)) return(wp)
  default_project
}

cdr             <- parse_bq_dataset(get_cdr_dataset())
billing_project <- get_billing_project(cdr$project)
message("Using AoU CDR dataset: ", cdr$path, " | billing project: ", billing_project)
fq_table <- function(table_name) paste0("`", cdr$path, ".", table_name, "`")
run_bq <- function(sql, bigint = "numeric") {
  job <- bq_project_query(x = billing_project, query = sql, use_legacy_sql = FALSE)
  bq_table_download(job, bigint = bigint, quiet = TRUE)
}
suppress_n <- function(d) d %>% mutate(across(starts_with("n_"), ~ ifelse(.x < SUPPRESS_MIN, NA_real_, .x)))

########################################################################
## 1. Structure of ds_survey
########################################################################

cols <- run_bq(paste0("SELECT column_name, data_type FROM `", cdr$path,
                      "`.INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'ds_survey' ORDER BY ordinal_position"))
write_excel_csv(cols, file.path(DATA_DIR, "T03_0_ds_survey_columns.csv"))
message("ds_survey columns: ", paste(cols$column_name, collapse = ", "))
has_col <- function(x) x %in% cols$column_name
date_col <- if (has_col("survey_datetime")) "survey_datetime" else NA_character_

########################################################################
## 2. Every survey: participants and date range
########################################################################

surveys <- run_bq(paste0(
  "SELECT survey, COUNT(DISTINCT person_id) AS n_persons",
  if (!is.na(date_col)) paste0(", MIN(DATE(", date_col, ")) AS first_date, MAX(DATE(", date_col, ")) AS last_date") else "",
  " FROM ", fq_table("ds_survey"), " GROUP BY survey ORDER BY n_persons DESC")) %>% suppress_n()
write_excel_csv(surveys, file.path(DATA_DIR, "T03_0_survey_names.csv"))
message("\n--- Surveys in ds_survey ---"); print(surveys, n = Inf)

hist_names <- surveys$survey[grepl(HISTORY_SURVEY_REGEX, surveys$survey, ignore.case = TRUE)]
if (length(hist_names) == 0)
  stop("No survey name matches '", HISTORY_SURVEY_REGEX, "'. Inspect T03_0_survey_names.csv and ",
       "set HISTORY_SURVEY_REGEX to the health-history survey name(s), then re-run.")
message("\nHealth-history survey(s): ", paste(hist_names, collapse = " | "))

########################################################################
## 3. Every question x answer of the health-history survey(s)
########################################################################

in_list <- paste0("'", gsub("'", "\\\\'", hist_names), "'", collapse = ", ")
version_cols <- intersect(c("survey_version_concept_id", "survey_version_name"), cols$column_name)
grp <- c("survey", version_cols, "question_concept_id", "question", "answer_concept_id", "answer")
items <- run_bq(paste0(
  "SELECT ", paste(grp, collapse = ", "), ", COUNT(DISTINCT person_id) AS n_persons",
  " FROM ", fq_table("ds_survey"),
  " WHERE survey IN (", in_list, ")",
  " GROUP BY ", paste(grp, collapse = ", "),
  " ORDER BY survey, question, answer")) %>% suppress_n()
write_excel_csv(items, file.path(DATA_DIR, "T03_0_health_history_items.csv"))
message("Question x answer rows: ", nrow(items), " | distinct questions: ", n_distinct(items$question_concept_id))

########################################################################
## 4. Items matching the crosswalk concepts
########################################################################

txt <- tolower(paste(items$question, items$answer))
cand <- bind_rows(lapply(names(CONCEPTS), function(k) {
  hit <- grepl(CONCEPTS[[k]], txt)
  if (!any(hit)) return(NULL)
  items[hit, ] %>% mutate(concept = k, .before = 1)
}))
missing_concepts <- setdiff(names(CONCEPTS), unique(cand$concept))
write_excel_csv(cand, file.path(DATA_DIR, "T03_0_crosswalk_candidates.csv"))

message("\n--- Crosswalk candidates: rows per concept ---")
print(cand %>% count(concept, name = "rows") %>% arrange(concept), n = Inf)
if (length(missing_concepts))
  message("No health-history item found for: ", paste(missing_concepts, collapse = ", "))
message("\nExamples (first 3 per concept):")
print(cand %>% group_by(concept) %>% slice_head(n = 3) %>% ungroup() %>%
        select(concept, question, answer, n_persons), n = Inf, width = Inf)

message("\nDone. Please send ", file.path(DATA_DIR, "T03_0_crosswalk_candidates.csv"),
        " (counts < ", SUPPRESS_MIN, " are already masked).")
