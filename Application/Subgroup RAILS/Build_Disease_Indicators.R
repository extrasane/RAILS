#### Build_Disease_Indicators.R — person x disease indicators (Workbench 2.0)
#### Run in the AoU Researcher Workbench (Verily / Workbench 2.0) workspace that
#### has the Controlled-Tier CDR attached as a resource. Uses the SAME CDR-
#### resolution grammar as 02_AoU_Prep.R (resource-resolved dataset, fully
#### qualified `project.dataset.table` references, run_bq) — NOT the old Terra
#### grammar (WORKSPACE_CDR / dbConnect to fc-aou-cdr-prod-ct), which is blocked
#### by VPC Service Controls after the Verily migration.
####
#### Case = a matching ICD9/10 code on >= 2 distinct dates (the "2 code
#### requirement"). Diseases with many ICD ranges are split into a ", pt. 2"
#### query and recombined, as in the reference notebook.
####
#### Only the GLOBAL person list + w_rails are attached here; the REGIONAL
#### sub-RAILS weight (w_subrails) is merged later in 10_State_Prevalence_Maps.R.
####
#### Before running, point the script at the CDR (one of):
####   Sys.setenv(AOU_CDR_RESOURCE_ID = "your-controlled-tier-cdr")   # Resources panel
####   Sys.setenv(AOU_CDR_DATASET     = "project.dataset")            # e.g. all-of-us-ehr-prod.C20XXQxRy
####
#### Inputs (../data/):  all_icd_codes.csv (ICD9CM/ICD10CM/Cause),
####                     global_rails_weights.csv (person_id, state, region, w_rails)
#### Output (../data/):  raking_wts_w_diseases_2_code_requirement.csv

library(tidyverse)
library(bigrquery)

########################################################################
## Workbench 2.0 CDR resolution (same grammar as 02_AoU_Prep.R)
########################################################################

AOU_CDR_RESOURCE_ID <- Sys.getenv("AOU_CDR_RESOURCE_ID", unset = "")
AOU_CDR_DATASET     <- Sys.getenv("AOU_CDR_DATASET",     unset = "")
BQ_BILLING_PROJECT  <- Sys.getenv("BQ_BILLING_PROJECT",  unset = "")

is_nonempty <- function(x) !is.na(x) && nzchar(trimws(x))

run_cmd <- function(cmd, args) {
  out <- tryCatch(system2(cmd, args = args, stdout = TRUE, stderr = TRUE),
                  error = function(e) character())
  out[!grepl("^$", out)]
}

resolve_wb_resource <- function(resource_id) {
  if (!is_nonempty(resource_id)) return(NA_character_)
  out <- run_cmd("wb", c("resolve", paste0("--id=", resource_id)))
  if (length(out) > 0 && any(grepl("\\.", out))) return(out[grepl("\\.", out)][1])
  out <- run_cmd("wb", c("resource", "resolve", paste0("--id=", resource_id)))
  if (length(out) > 0 && any(grepl("\\.", out))) return(out[grepl("\\.", out)][1])
  stop("Could not resolve Workbench resource id: ", resource_id,
       "\nCheck that the resource is attached to this workspace and `wb` is available.")
}

get_cdr_dataset <- function() {
  if (is_nonempty(AOU_CDR_DATASET)) return(AOU_CDR_DATASET)
  legacy <- Sys.getenv("WORKSPACE_CDR", unset = "")
  if (is_nonempty(legacy)) return(legacy)
  if (is_nonempty(AOU_CDR_RESOURCE_ID)) return(resolve_wb_resource(AOU_CDR_RESOURCE_ID))
  stop("No AoU CDR dataset found. Set AOU_CDR_RESOURCE_ID or AOU_CDR_DATASET ",
       "(see the Workbench 2.0 Resources panel).")
}

parse_bq_dataset <- function(dataset_path) {
  dataset_path <- gsub("`", "", dataset_path)
  parts <- strsplit(dataset_path, "\\.")[[1]]
  if (length(parts) != 2) stop("Expected 'project.dataset', got: ", dataset_path)
  list(project = parts[1], dataset = parts[2], path = paste(parts[1], parts[2], sep = "."))
}

get_billing_project <- function(default_project) {
  candidates <- c(BQ_BILLING_PROJECT,
                  Sys.getenv("GOOGLE_PROJECT",       unset = ""),
                  Sys.getenv("GOOGLE_CLOUD_PROJECT", unset = ""),
                  Sys.getenv("GCLOUD_PROJECT",       unset = ""))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) > 0) return(candidates[1])
  gp <- run_cmd("gcloud", c("config", "get-value", "project"))
  gp <- gp[nzchar(gp) & !grepl("^Your active configuration", gp)]
  if (length(gp) > 0 && !grepl("ERROR|unset", gp[1], ignore.case = TRUE)) return(gp[1])
  default_project
}

cdr             <- parse_bq_dataset(get_cdr_dataset())
billing_project <- get_billing_project(cdr$project)
message("Using AoU CDR dataset: ", cdr$path, " | billing: ", billing_project)

fq_table <- function(table_name) paste0("`", cdr$path, ".", table_name, "`")

## person_id is stored as double throughout this pipeline, so download BIGINT
## as numeric to keep the %in% join consistent.
run_bq <- function(sql, bigint = "numeric") {
  job <- bq_project_query(x = billing_project, query = sql, use_legacy_sql = FALSE)
  bq_table_download(job, bigint = bigint, quiet = TRUE)
}

########################################################################
## Person base: global RAILS weights (../data/)
########################################################################

raking.wts <- read_csv("../data/global_rails_weights.csv", show_col_types = FALSE) %>%
  select(person_id, state, region, w_rails)
message("Global-weighted participants: ", nrow(raking.wts))

########################################################################
## Disease -> ICD lookup (drop blank rows, rename to icd9 / icd10 / Disease)
########################################################################

codes <- read_csv("../data/all_icd_codes.csv", show_col_types = FALSE) %>%
  filter(!(is.na(ICD10CM) & is.na(ICD9CM))) %>%
  rename(icd9 = ICD9CM, icd10 = ICD10CM, Disease = Cause)
message("Diseases in lookup: ", nrow(codes))

########################################################################
## Helpers (from the reference notebook)
########################################################################

reorg <- function(vec) {
  tibble(code = vec) %>%
    tidyr::separate(code, into = c("first", "last"), sep = "-", fill = "right") %>%
    mutate(last = ifelse(is.na(last), first, last)) %>%
    distinct()
}

## 2-code GROUP BY Person query for one disease, with FULLY QUALIFIED tables
text.creator <- function(df) {
  clauses <- sprintf("(condition_source_value BETWEEN '%s' AND '%s')",
                     df$first, df$last)
  txt <- paste(clauses, collapse = " OR ")
  str_glue("
    SELECT person_id AS Person,
           COUNT(DISTINCT condition_start_date) AS N
    FROM {fq_table('condition_occurrence')} AS co
    JOIN {fq_table('concept')} AS mc ON (condition_source_concept_id = mc.concept_id)
    WHERE mc.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    AND ({txt})
    GROUP BY Person")
}

########################################################################
## Assemble per-disease code lists; split diseases with > 100 ranges
########################################################################

icd9  <- lapply(strsplit(codes$icd9,  split = ", "), reorg)
icd10 <- lapply(strsplit(codes$icd10, split = ", "), reorg)

codes.to.search <- lapply(seq_len(nrow(codes)),
                          function(i) na.omit(rbind(icd9[[i]], icd10[[i]])))
names(codes.to.search) <- codes$Disease

n <- length(codes.to.search)
for (j in seq_len(n)) {
  temp <- codes.to.search[[j]]
  if (nrow(temp) > 100) {
    N <- length(codes.to.search)
    codes.to.search[[N + 1]]      <- temp[76:nrow(temp), ]
    codes.to.search[[j]]          <- temp[1:75, ]
    names(codes.to.search)[N + 1] <- paste0(names(codes.to.search)[j], ", pt. 2")
  }
}

########################################################################
## Query each disease (run_bq), keep 2-code cases, recombine split parts
########################################################################

results       <- lapply(lapply(codes.to.search, text.creator), run_bq)
names(results) <- names(codes.to.search)
results.2code <- lapply(results, function(x) filter(x, N > 1))

for (nm2 in grep(", pt\\. 2$", names(results.2code), value = TRUE)) {
  base <- sub(", pt\\. 2$", "", nm2)
  results.2code[[base]] <- rbind(results.2code[[base]], results.2code[[nm2]])
}

########################################################################
## Assign a 0/1 indicator column per disease, then save
########################################################################

for (i in seq_len(nrow(codes))) {
  dis <- codes$Disease[i]
  raking.wts[[dis]] <- as.integer(raking.wts$person_id %in% results.2code[[dis]]$Person)
  message("  ", dis, ": ", sum(raking.wts[[dis]]), " cases")
}

write_excel_csv(raking.wts, "../data/raking_wts_w_diseases_2_code_requirement.csv")
