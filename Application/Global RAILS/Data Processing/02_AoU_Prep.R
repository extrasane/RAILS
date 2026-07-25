#### 02_AoU_Prep.R — Query and prepare AoU biobank data
#### Run inside the AoU Researcher Workbench
####
#### The Workbench 1.0 implementation is retained, commented out, below the
#### header for reference. The active code uses the Workbench 2.0 grammar:
#### resource-resolved CDR, direct bq_project_query downloads, and fully
#### qualified `project.dataset.table` references.

library(tidyverse)
library(bigrquery)
library(haven)
library(httr)
library(jsonlite)
library(tidyr)
library(lubridate)
library(survey)

# ########################################################################
# ## WORKBENCH 1.0 (retained for reference — not run)
# ########################################################################
#
# ########################################################################
# ## BigQuery queries
# ## Replace YOUR_DATASET with your workspace's controlled-tier dataset
# ## identifier (visible in the Workbench dataset panel under "Controlled CDR")
# ########################################################################
#
# ## All queries restrict the cohort to participants with EHR data
# ## (has_ehr_data = 1), matching the "Raking Weighting v1" dataset definition.
# ehr_cohort_sql <- "
#     SELECT distinct person_id
#     FROM `cb_search_person` cb_search_person
#     WHERE cb_search_person.person_id IN (
#       SELECT person_id FROM `cb_search_person` p WHERE has_ehr_data = 1
#     )"
#
# ## --- Demographics ---
# dataset_person_sql <- paste0("
#   SELECT
#     person.person_id,
#     p_gender_concept.concept_name    AS gender,
#     person.birth_datetime            AS date_of_birth,
#     p_race_concept.concept_name      AS race,
#     p_ethnicity_concept.concept_name AS ethnicity,
#     p_sex_at_birth_concept.concept_name AS sex_at_birth
#   FROM `person` person
#   LEFT JOIN `concept` p_gender_concept
#     ON person.gender_concept_id = p_gender_concept.concept_id
#   LEFT JOIN `concept` p_race_concept
#     ON person.race_concept_id = p_race_concept.concept_id
#   LEFT JOIN `concept` p_ethnicity_concept
#     ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id
#   LEFT JOIN `concept` p_sex_at_birth_concept
#     ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
#   WHERE person.PERSON_ID IN (", ehr_cohort_sql, ")")
#
# ## --- Survey responses (income, homeown, edu, care place) ---
# dataset_survey_sql <- paste0("
#   SELECT
#     answer.person_id,
#     answer.question_concept_id,
#     answer.question,
#     answer.answer
#   FROM `ds_survey` answer
#   WHERE question_concept_id IN (1585370, 1585375, 1585899, 1585940, 43530593)
#     AND answer.PERSON_ID IN (", ehr_cohort_sql, ")")
#
# ## --- State of residence ---
# dataset_state_sql <- paste0("
#   SELECT
#     person.person_id,
#     person.state_of_residence_source_value AS state
#   FROM `person_ext` person
#   WHERE person.PERSON_ID IN (", ehr_cohort_sql, ")")
#
# ########################################################################
# ## Export each query from BigQuery to Cloud Storage, then read back
# ## NOTE: bq_table_save only needs to run once per day.
# ##       Re-running on the same day overwrites; a new day writes a new path.
# ########################################################################
#
# export_bq <- function(sql, export_name, col_types) {
#   export_path <- file.path(
#     Sys.getenv("WORKSPACE_BUCKET"), "bq_exports",
#     Sys.getenv("OWNER_EMAIL"),
#     strftime(lubridate::now(), "%Y%m%d"),
#     export_name,
#     paste0(export_name, "_*.csv")
#   )
#   message(str_glue("Exporting {export_name} to {export_path}"))
#   bq_table_save(
#     bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), sql, billing = Sys.getenv("GOOGLE_PROJECT")),
#     export_path,
#     destination_format = "CSV"
#   )
#   bind_rows(
#     map(system2("gsutil", args = c("ls", export_path), stdout = TRUE, stderr = TRUE),
#         function(csv) {
#           message(str_glue("Loading {csv}."))
#           chunk <- read_csv(pipe(str_glue("gsutil cat {csv}")),
#                             col_types = col_types, show_col_types = FALSE)
#           chunk
#         })
#   )
# }
#
# person_df <- export_bq(
#   dataset_person_sql, "person",
#   cols(gender = col_character(), race = col_character(),
#        ethnicity = col_character(), sex_at_birth = col_character())
# )
#
# survey_df <- export_bq(
#   dataset_survey_sql, "survey",
#   cols(question = col_character(), answer = col_character())
# )
#
# state_df <- export_bq(
#   dataset_state_sql, "person_state",
#   cols(state = col_character())
# )
#
# ########################################################################
# ## Join tables and recode variables
# ########################################################################
#
# ## values_fn = first guards against participants with multiple answers to the
# ## same question (e.g. re-taken surveys), which would otherwise produce
# ## list-columns and break the downstream recoding.
# survey_wide <- survey_df %>%
#   select(-question_concept_id) %>%
#   filter(question != "The Basics: Sexual Orientation") %>%
#   pivot_wider(names_from = question, values_from = answer, values_fn = first)
#
# dt0 <- person_df %>%
#   full_join(survey_wide, by = "person_id") %>%
#   full_join(state_df,    by = "person_id") %>%
#   rename(
#     income    = `Income: Annual Income`,
#     homeown   = `Home Own: Current Home Own`,
#     edu       = `Education Level: Highest Grade`,
#     careplace = `Health Advice: What Kind Of Place`
#   )
#
# dt0 <- dt0 %>%
#   mutate(date_of_birth = as.Date(date_of_birth)) %>%
#   mutate(age = round(interval(date_of_birth, as.Date("2024-08-01")) /
#                        duration(num = 1, units = "years"), digits = 0)) %>%
#   select(-date_of_birth) %>%
#   mutate(race = case_when(
#     race %in% c("I prefer not to answer", "PMI: Skip")                        ~ NA_character_,
#     race %in% c("None Indicated", "None of these",
#                 "Middle Eastern or North African",
#                 "Native Hawaiian or Other Pacific Islander",
#                 "More than one population")                                    ~ "Others",
#     TRUE                                                                       ~ race
#   )) %>%
#   mutate(ethnicity = case_when(
#     ethnicity == "Hispanic or Latino"                                          ~ "Yes",
#     ethnicity %in% c("No matching concept",
#                      "Not Hispanic or Latino",
#                      "What Race Ethnicity: Race Ethnicity None Of These")      ~ "No",
#     ethnicity %in% c("PMI: Prefer Not To Answer", "PMI: Skip")                ~ NA_character_
#   )) %>%
#   mutate(sex = case_when(
#     sex_at_birth %in% c("I prefer not to answer", "Intersex",
#                         "No matching concept", "None", "PMI: Skip")            ~ NA_character_,
#     TRUE                                                                       ~ sex_at_birth
#   )) %>%
#   select(-sex_at_birth) %>%
#   mutate(income = case_when(
#     income %in% c("Annual Income: 100k 150k",
#                   "Annual Income: 150k 200k",
#                   "Annual Income: more 200k")                                  ~ ">100k",
#     income == "Annual Income: 75k 100k"                                        ~ "75k-100k",
#     income == "Annual Income: 50k 75k"                                         ~ "50k-75k",
#     income == "Annual Income: 35k 50k"                                         ~ "35k-50k",
#     income %in% c("Annual Income: 25k 35k",
#                   "Annual Income: 10k 25k",
#                   "Annual Income: less 10k")                                   ~ "<35k",
#     TRUE                                                                       ~ NA_character_
#   )) %>%
#   ## "Other" (not "Others") here; the
#   ## analysis step recodes it to "Others" via ifelse before setting levels.
#   mutate(homeown = case_when(
#     homeown %in% c("Current Home Own: Other Arrangement",
#                    "PMI: Dont Know")                                           ~ "Other",
#     homeown == "Current Home Own: Own"                                         ~ "Own",
#     homeown == "Current Home Own: Rent"                                        ~ "Rent",
#     TRUE                                                                       ~ NA_character_
#   )) %>%
#   mutate(edu = case_when(
#     edu %in% c("Highest Grade: Advanced Degree",
#                "Highest Grade: College Graduate")                              ~ "College graduate or advanced",
#     edu == "Highest Grade: College One to Three"                               ~ "Some college",
#     edu == "Highest Grade: Twelve Or GED"                                      ~ "Highschool graduate",
#     edu == "Highest Grade: Nine Through Eleven"                                ~ "Some highschool",
#     edu %in% c("Highest Grade: Never Attended",
#                "Highest Grade: One Through Four",
#                "Highest Grade: Five Through Eight")                            ~ "Less than highschool",
#     TRUE                                                                       ~ NA_character_
#   )) %>%
#   mutate(careplace = case_when(
#     careplace == "What Kind Of Place: Doctors Office"          ~ "Doctors office",
#     careplace == "What Kind Of Place: Emergency Room"          ~ "Emergency",
#     careplace == "What Kind Of Place: No One Place Most Often" ~ "None",
#     careplace == "What Kind Of Place: Urgent Care"             ~ "Urgent",
#     careplace == "What Kind Of Place: Some Other Place"        ~ "Others",
#     TRUE                                                       ~ NA_character_
#   ))
#
# ## Extract two-letter state abbreviation from "PII State: XX" format
# dt00 <- dt0 %>%
#   mutate(state = ifelse(str_detect(state, "^PII State: "),
#                         str_extract(state, "(?<=PII State: )\\w{2}"),
#                         NA))
#
# ########################################################################
# ## Save raw cleaned AoU data to bucket
# ########################################################################
#
# write_excel_csv(dt00, "aou_raking_dt.csv")
# my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
# system(paste0("gsutil cp ./aou_raking_dt.csv ", my_bucket, "/data/"), intern = TRUE)
#
# ########################################################################
# ## Prepare individual-level and aggregated AoU datasets
# ########################################################################
#
# ## Define the states in each region
# south     <- c("AL","AR","FL","GA","KY","LA","MS","NC","SC","TN","TX","VA","WV")
# midwest   <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI")
# northeast <- c("CT","DE","ME","MD","MA","NH","NJ","NY","PA","RI","VT")
# west      <- c("AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY")
#
# ## Create a named vector to map states to regions
# state_region_map <- c(
#   setNames(rep("South",     length(south)),     south),
#   setNames(rep("Midwest",   length(midwest)),   midwest),
#   setNames(rep("Northeast", length(northeast)), northeast),
#   setNames(rep("West",      length(west)),      west)
# )
#
# dt_aou <- dt00 %>%
#   select(!c(gender, careplace)) %>%
#   na.omit()
#
# dt_aou <- dt_aou %>%
#   mutate(race_eth = case_when(
#     ethnicity == "Yes"                       ~ "Hispanic",
#     race      == "Asian"                     ~ "NH Asian",
#     race      == "Black or African American" ~ "NH Black",
#     race      == "White"                     ~ "NH White",
#     TRUE                                     ~ "Others"
#   ))
#
# dt_aou <- dt_aou %>%
#   mutate(agegroup = case_when(
#     age <= 24            ~ "18-24",
#     age > 24 & age <= 44 ~ "25-44",
#     age > 44 & age <= 64 ~ "45-64",
#     age > 64 & age <= 74 ~ "65-74",
#     TRUE                 ~ "75+"
#   ))
#
# dt_aou <- dt_aou %>%
#   mutate(sex      = factor(sex,      levels = c("Female","Male"))) %>%
#   mutate(race_eth = factor(race_eth, levels = c("Hispanic","NH Asian","NH Black","NH White","Others"))) %>%
#   mutate(income   = factor(income,   levels = c("<35k","35k-50k","50k-75k","75k-100k",">100k"))) %>%
#   mutate(agegroup = factor(agegroup, levels = c("18-24","25-44","45-64","65-74","75+"))) %>%
#   mutate(edu      = factor(edu,      levels = c("Less than highschool","Some highschool",
#                                                 "Highschool graduate","Some college",
#                                                 "College graduate or advanced"))) %>%
#   mutate(homeown  = ifelse(homeown == "Other", "Others", homeown)) %>%
#   mutate(homeown  = factor(homeown,  levels = c("Own","Rent","Others"))) %>%
#   na.omit()
#
# dt_aou <- dt_aou %>%
#   mutate(region = state_region_map[state]) %>%
#   mutate(region = factor(region, levels = c("Northeast","Midwest","South","West"))) %>%
#   na.omit() %>%
#   mutate(weight = 1)
#
# names_univar <- c("agegroup","sex","edu","homeown","income","race_eth","region")
# cat_formula  <- formula(paste0("weight ~", paste(names_univar, collapse = "+")))
# dt_agg_aou   <- aggregate(cat_formula, data = dt_aou, sum)
#
# write_excel_csv(dt_agg_aou, "dt_agg_aou_v3.csv")
# system(paste0("gsutil cp ./dt_agg_aou_v3.csv ", my_bucket, "/data/"), intern = TRUE)
# system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = TRUE)


########################################################################
## WORKBENCH 2.0 (active)
########################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(bigrquery)
  library(lubridate)
})

# ========================================================================
# User configuration
# ========================================================================

# Option A, recommended in Workbench 2.0:
# In the R terminal or notebook, set this to the resource id/name shown in Resources.
# Example:
#   Sys.setenv(AOU_CDR_RESOURCE_ID = "my-controlled-tier-cdr")
AOU_CDR_RESOURCE_ID <- Sys.getenv("AOU_CDR_RESOURCE_ID", unset = "")

# Option B, direct BigQuery dataset path:
# Example:
#   Sys.setenv(AOU_CDR_DATASET = "all-of-us-ehr-prod.C2025Q4R6")
AOU_CDR_DATASET <- Sys.getenv("AOU_CDR_DATASET", unset = "")

# Optional: billing project. If blank, the script tries common env vars, then
# falls back to the project part of the CDR dataset.
BQ_BILLING_PROJECT <- Sys.getenv("BQ_BILLING_PROJECT", unset = "")

# Optional: copy final CSV files to a GCS bucket resource.
# Example:
#   Sys.setenv(OUTPUT_BUCKET_RESOURCE_ID = "my-output-bucket")
OUTPUT_BUCKET_RESOURCE_ID <- Sys.getenv("OUTPUT_BUCKET_RESOURCE_ID", unset = "")
OUTPUT_BUCKET_PATH        <- Sys.getenv("OUTPUT_BUCKET_PATH", unset = "data")

# Age reference date. Keep fixed so results are reproducible across runs.
AGE_REFERENCE_DATE <- as.Date(Sys.getenv("AGE_REFERENCE_DATE", unset = "2024-08-01"))

# ========================================================================
# Helpers for Workbench 2.0 / BigQuery
# ========================================================================

is_nonempty <- function(x) !is.na(x) && nzchar(trimws(x))

run_cmd <- function(cmd, args) {
  out <- tryCatch(
    system2(cmd, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  out <- out[!grepl("^$", out)]
  out
}

resolve_wb_resource <- function(resource_id) {
  if (!is_nonempty(resource_id)) return(NA_character_)

  # Newer Workbench CLI supports both `wb resolve --id=...` and
  # `wb resource resolve --id=...`; try both to be robust across images.
  out <- run_cmd("wb", c("resolve", paste0("--id=", resource_id)))
  if (length(out) > 0 && any(grepl("\\.", out))) return(out[grepl("\\.", out)][1])

  out <- run_cmd("wb", c("resource", "resolve", paste0("--id=", resource_id)))
  if (length(out) > 0 && any(grepl("\\.", out))) return(out[grepl("\\.", out)][1])

  stop(
    "Could not resolve Workbench resource id: ", resource_id, "\n",
    "Check that the resource is attached to this workspace and that `wb` is available."
  )
}

get_cdr_dataset <- function() {
  if (is_nonempty(AOU_CDR_DATASET)) return(AOU_CDR_DATASET)

  # Legacy fallback. Useful only if AoU still populates it in your image.
  legacy <- Sys.getenv("WORKSPACE_CDR", unset = "")
  if (is_nonempty(legacy)) return(legacy)

  if (is_nonempty(AOU_CDR_RESOURCE_ID)) return(resolve_wb_resource(AOU_CDR_RESOURCE_ID))

  stop(
    "No AoU CDR dataset found.\n\n",
    "Set one of the following before running this script:\n",
    "  Sys.setenv(AOU_CDR_RESOURCE_ID = 'your-data-collection-resource-id')\n",
    "or\n",
    "  Sys.setenv(AOU_CDR_DATASET = 'project.dataset')\n\n",
    "You can usually find the resource id/name in the Workbench 2.0 Resources panel."
  )
}

parse_bq_dataset <- function(dataset_path) {
  dataset_path <- gsub("`", "", dataset_path)
  parts <- strsplit(dataset_path, "\\.")[[1]]
  if (length(parts) != 2) {
    stop("Expected BigQuery dataset path as 'project.dataset', got: ", dataset_path)
  }
  list(project = parts[1], dataset = parts[2], path = paste(parts[1], parts[2], sep = "."))
}

get_billing_project <- function(default_project) {
  candidates <- c(
    BQ_BILLING_PROJECT,
    Sys.getenv("GOOGLE_PROJECT", unset = ""),
    Sys.getenv("GOOGLE_CLOUD_PROJECT", unset = ""),
    Sys.getenv("GCLOUD_PROJECT", unset = "")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) > 0) return(candidates[1])

  gcloud_project <- run_cmd("gcloud", c("config", "get-value", "project"))
  gcloud_project <- gcloud_project[!grepl("^Your active configuration", gcloud_project)]
  gcloud_project <- gcloud_project[nzchar(gcloud_project)]
  if (length(gcloud_project) > 0 && !grepl("ERROR|unset", gcloud_project[1], ignore.case = TRUE)) {
    return(gcloud_project[1])
  }

  default_project
}

cdr <- parse_bq_dataset(get_cdr_dataset())
billing_project <- get_billing_project(cdr$project)

message("Using AoU CDR dataset: ", cdr$path)
message("Using BigQuery billing project: ", billing_project)

fq_table <- function(table_name) {
  paste0("`", cdr$path, ".", table_name, "`")
}

table_exists <- function(table_name) {
  ok <- tryCatch({
    bq_table_exists(bq_table(cdr$project, cdr$dataset, table_name))
  }, error = function(e) FALSE)
  isTRUE(ok)
}

run_bq <- function(sql, bigint = "integer64") {
  message("Submitting BigQuery job...")
  job <- bq_project_query(
    x = billing_project,
    query = sql,
    use_legacy_sql = FALSE
  )
  bq_table_download(job, bigint = bigint, quiet = TRUE)
}

# ========================================================================
# Cohort definition
# ========================================================================

# In legacy AoU/Terra notebooks, `cb_search_person` was often available and
# contained `has_ehr_data`. In Workbench 2.0/Data Explorer, this table may or
# may not exist depending on how the Data Collection was attached/generated.
# This script uses it when present; otherwise it falls back to all persons and
# prints a clear warning. If your new workspace has a generated cohort table,
# replace ehr_cohort_sql with that table.

if (table_exists("cb_search_person")) {
  message("Found cb_search_person; restricting to has_ehr_data = 1.")
  ehr_cohort_sql <- paste0(
    "SELECT DISTINCT person_id\n",
    "FROM ", fq_table("cb_search_person"), "\n",
    "WHERE has_ehr_data = 1"
  )
} else {
  warning(
    "Table cb_search_person was not found in ", cdr$path, ".\n",
    "Proceeding with all rows in person. If you need the exact legacy EHR-only cohort, ",
    "create/export the cohort in Data Explorer or identify the new cohort table and ",
    "replace ehr_cohort_sql manually."
  )
  ehr_cohort_sql <- paste0(
    "SELECT DISTINCT person_id\n",
    "FROM ", fq_table("person")
  )
}

# ========================================================================
# BigQuery queries
# ========================================================================

# --- Demographics ---
dataset_person_sql <- paste0(
"SELECT
  person.person_id,
  p_gender_concept.concept_name AS gender,
  person.birth_datetime AS date_of_birth,
  p_race_concept.concept_name AS race,
  p_ethnicity_concept.concept_name AS ethnicity,
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
WHERE person.person_id IN (", ehr_cohort_sql, ")"
)

# --- Survey responses: income, homeown, education, care place, sexual orientation excluded later ---
if (!table_exists("ds_survey")) {
  stop(
    "Table ds_survey was not found in ", cdr$path, ".\n",
    "In Workbench 2.0, you may need to add the AoU survey Data Collection/resource ",
    "or use the survey table generated by Data Explorer."
  )
}

dataset_survey_sql <- paste0(
"SELECT
  answer.person_id,
  answer.question_concept_id,
  answer.question,
  answer.answer
FROM ", fq_table("ds_survey"), " answer
WHERE question_concept_id IN (1585370, 1585375, 1585899, 1585940, 43530593)
  AND answer.person_id IN (", ehr_cohort_sql, ")"
)

# --- State of residence ---
if (!table_exists("person_ext")) {
  stop(
    "Table person_ext was not found in ", cdr$path, ".\n",
    "Your state extraction depends on person_ext.state_of_residence_source_value. ",
    "Check whether the Controlled Tier Data Collection is attached."
  )
}

dataset_state_sql <- paste0(
"SELECT
  person.person_id,
  person.state_of_residence_source_value AS state
FROM ", fq_table("person_ext"), " person
WHERE person.person_id IN (", ehr_cohort_sql, ")"
)

# ========================================================================
# Execute queries directly into R
# ========================================================================

person_df <- run_bq(dataset_person_sql) %>%
  mutate(
    gender = as.character(gender),
    race = as.character(race),
    ethnicity = as.character(ethnicity),
    sex_at_birth = as.character(sex_at_birth)
  )

survey_df <- run_bq(dataset_survey_sql) %>%
  mutate(
    question = as.character(question),
    answer = as.character(answer)
  )

state_df <- run_bq(dataset_state_sql) %>%
  mutate(state = as.character(state))

message("Rows downloaded: person=", nrow(person_df),
        ", survey=", nrow(survey_df),
        ", state=", nrow(state_df))

# ========================================================================
# Join tables and recode variables
# ========================================================================

survey_wide <- survey_df %>%
  select(-question_concept_id) %>%
  filter(question != "The Basics: Sexual Orientation") %>%
  pivot_wider(names_from = question, values_from = answer, values_fn = first)

required_survey_cols <- c(
  "Income: Annual Income",
  "Home Own: Current Home Own",
  "Education Level: Highest Grade",
  "Health Advice: What Kind Of Place"
)
missing_cols <- setdiff(required_survey_cols, names(survey_wide))
if (length(missing_cols) > 0) {
  stop(
    "The following expected survey columns were not found after pivot_wider():\n  ",
    paste(missing_cols, collapse = "\n  "),
    "\nCheck question_concept_id values and answer/question labels in the new CDR."
  )
}

dt0 <- person_df %>%
  full_join(survey_wide, by = "person_id") %>%
  full_join(state_df, by = "person_id") %>%
  rename(
    income    = `Income: Annual Income`,
    homeown   = `Home Own: Current Home Own`,
    edu       = `Education Level: Highest Grade`,
    careplace = `Health Advice: What Kind Of Place`
  )

dt0 <- dt0 %>%
  mutate(date_of_birth = as.Date(date_of_birth)) %>%
  mutate(age = round(interval(date_of_birth, AGE_REFERENCE_DATE) /
                       duration(num = 1, units = "years"), digits = 0)) %>%
  select(-date_of_birth) %>%
  mutate(race = case_when(
    race %in% c("I prefer not to answer", "PMI: Skip") ~ NA_character_,
    race %in% c("None Indicated", "None of these",
                "Middle Eastern or North African",
                "Native Hawaiian or Other Pacific Islander",
                "More than one population") ~ "Others",
    TRUE ~ race
  )) %>%
  mutate(ethnicity = case_when(
    ethnicity == "Hispanic or Latino" ~ "Yes",
    ethnicity %in% c("No matching concept",
                     "Not Hispanic or Latino",
                     "What Race Ethnicity: Race Ethnicity None Of These") ~ "No",
    ethnicity %in% c("PMI: Prefer Not To Answer", "PMI: Skip") ~ NA_character_,
    TRUE ~ ethnicity
  )) %>%
  mutate(sex = case_when(
    sex_at_birth %in% c("I prefer not to answer", "Intersex",
                        "No matching concept", "None", "PMI: Skip") ~ NA_character_,
    TRUE ~ sex_at_birth
  )) %>%
  select(-sex_at_birth) %>%
  mutate(income = case_when(
    income %in% c("Annual Income: 100k 150k",
                  "Annual Income: 150k 200k",
                  "Annual Income: more 200k") ~ ">100k",
    income == "Annual Income: 75k 100k" ~ "75k-100k",
    income == "Annual Income: 50k 75k" ~ "50k-75k",
    income == "Annual Income: 35k 50k" ~ "35k-50k",
    income %in% c("Annual Income: 25k 35k",
                  "Annual Income: 10k 25k",
                  "Annual Income: less 10k") ~ "<35k",
    TRUE ~ NA_character_
  )) %>%
  mutate(homeown = case_when(
    homeown %in% c("Current Home Own: Other Arrangement",
                   "PMI: Dont Know") ~ "Other",
    homeown == "Current Home Own: Own" ~ "Own",
    homeown == "Current Home Own: Rent" ~ "Rent",
    TRUE ~ NA_character_
  )) %>%
  mutate(edu = case_when(
    edu %in% c("Highest Grade: Advanced Degree",
               "Highest Grade: College Graduate") ~ "College graduate or advanced",
    edu == "Highest Grade: College One to Three" ~ "Some college",
    edu == "Highest Grade: Twelve Or GED" ~ "Highschool graduate",
    edu == "Highest Grade: Nine Through Eleven" ~ "Some highschool",
    edu %in% c("Highest Grade: Never Attended",
               "Highest Grade: One Through Four",
               "Highest Grade: Five Through Eight") ~ "Less than highschool",
    TRUE ~ NA_character_
  )) %>%
  mutate(careplace = case_when(
    careplace == "What Kind Of Place: Doctors Office" ~ "Doctors office",
    careplace == "What Kind Of Place: Emergency Room" ~ "Emergency",
    careplace == "What Kind Of Place: No One Place Most Often" ~ "None",
    careplace == "What Kind Of Place: Urgent Care" ~ "Urgent",
    careplace == "What Kind Of Place: Some Other Place" ~ "Others",
    TRUE ~ NA_character_
  ))

# Extract two-letter state abbreviation from "PII State: XX" format.
dt00 <- dt0 %>%
  mutate(state = ifelse(str_detect(state, "^PII State: "),
                        str_extract(state, "(?<=PII State: )\\w{2}"),
                        NA_character_))

# ========================================================================
# Prepare individual-level and aggregated AoU datasets
# ========================================================================

south     <- c("AL", "AR", "FL", "GA", "KY", "LA", "MS", "NC", "SC", "TN", "TX", "VA", "WV")
midwest   <- c("IL", "IN", "IA", "KS", "MI", "MN", "MO", "NE", "ND", "OH", "SD", "WI")
northeast <- c("CT", "DE", "ME", "MD", "MA", "NH", "NJ", "NY", "PA", "RI", "VT")
west      <- c("AK", "AZ", "CA", "CO", "HI", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")

state_region_map <- c(
  setNames(rep("South", length(south)), south),
  setNames(rep("Midwest", length(midwest)), midwest),
  setNames(rep("Northeast", length(northeast)), northeast),
  setNames(rep("West", length(west)), west)
)

dt_aou <- dt00 %>%
  select(!c(gender, careplace)) %>%
  na.omit() %>%
  mutate(race_eth = case_when(
    ethnicity == "Yes" ~ "Hispanic",
    race == "Asian" ~ "NH Asian",
    race == "Black or African American" ~ "NH Black",
    race == "White" ~ "NH White",
    TRUE ~ "Others"
  )) %>%
  mutate(agegroup = case_when(
    age <= 24 ~ "18-24",
    age > 24 & age <= 44 ~ "25-44",
    age > 44 & age <= 64 ~ "45-64",
    age > 64 & age <= 74 ~ "65-74",
    TRUE ~ "75+"
  )) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male"))) %>%
  mutate(race_eth = factor(race_eth, levels = c("Hispanic", "NH Asian", "NH Black", "NH White", "Others"))) %>%
  mutate(income = factor(income, levels = c("<35k", "35k-50k", "50k-75k", "75k-100k", ">100k"))) %>%
  mutate(agegroup = factor(agegroup, levels = c("18-24", "25-44", "45-64", "65-74", "75+"))) %>%
  mutate(edu = factor(edu, levels = c("Less than highschool", "Some highschool",
                                      "Highschool graduate", "Some college",
                                      "College graduate or advanced"))) %>%
  mutate(homeown = ifelse(homeown == "Other", "Others", homeown)) %>%
  mutate(homeown = factor(homeown, levels = c("Own", "Rent", "Others"))) %>%
  na.omit() %>%
  mutate(region = state_region_map[state]) %>%
  mutate(region = factor(region, levels = c("Northeast", "Midwest", "South", "West"))) %>%
  na.omit() %>%
  mutate(weight = 1)

names_univar <- c("agegroup", "sex", "edu", "homeown", "income", "race_eth", "region")
cat_formula  <- formula(paste0("weight ~ ", paste(names_univar, collapse = "+")))
dt_agg_aou   <- aggregate(cat_formula, data = dt_aou, sum)

# ========================================================================
# Save outputs locally and optionally copy to GCS
# ========================================================================

write_excel_csv(dt00, "aou_raking_dt.csv")
write_excel_csv(dt_aou, "dt_aou_individual_v3.csv")
write_excel_csv(dt_agg_aou, "dt_agg_aou_v3.csv")

copy_to_gcs_if_requested <- function(files) {
  if (!is_nonempty(OUTPUT_BUCKET_RESOURCE_ID)) {
    message("OUTPUT_BUCKET_RESOURCE_ID not set; files saved locally only.")
    return(invisible(NULL))
  }

  bucket <- resolve_wb_resource(OUTPUT_BUCKET_RESOURCE_ID)
  bucket <- sub("/$", "", bucket)
  dest <- paste0(bucket, "/", OUTPUT_BUCKET_PATH, "/")
  message("Copying outputs to ", dest)

  for (f in files) {
    status <- system2("gsutil", args = c("cp", f, dest), stdout = TRUE, stderr = TRUE)
    message(paste(status, collapse = "\n"))
  }

  invisible(NULL)
}

copy_to_gcs_if_requested(c("aou_raking_dt.csv", "dt_aou_individual_v3.csv", "dt_agg_aou_v3.csv"))

message("Done.")
message("Local outputs:")
message("  - aou_raking_dt.csv")
message("  - dt_aou_individual_v3.csv")
message("  - dt_agg_aou_v3.csv")
