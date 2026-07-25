#### 09_Sub_RAILS.R — Single-subgroup RAILS analysis
#### Run inside the AoU Researcher Workbench.
#### Requires AoU_Fun.R (from ../Global RAILS/RAILS Procedure/) in the same
#### directory.
####
#### REUSABLE for any single stratum: set SUB_VAR / SUB_LEVEL below (e.g.
#### sex == "Female" for the breast cancer study, or region == "South").
#### The subgroup variable is dropped from the model since it is constant
#### within the stratum, and weights are scaled to that subgroup's PUMS
#### total. Because only one stratum is analysed, this calls
#### fun.rails.threeway directly rather than the by-level wrapper
#### fun.sub.rails.threeway (use that one to run every level at once).

library(tidyverse)
library(survey)
library(Matrix)
source("AoU_Fun.R")

########################################################################
## SUBGROUP SETTING — the only lines to change between analyses
########################################################################

SUB_VAR   <- "sex"        # stratification variable
SUB_LEVEL <- "Female"     # level to analyse
OUT_FILE  <- "dt_sub_aou_femaleonly.csv"

## All seven shared covariates; SUB_VAR is removed automatically below
names_univar_all <- c("agegroup", "edu", "homeown", "income",
                      "race_eth", "sex", "region")

########################################################################
## Load datasets from Workspace bucket
########################################################################

my_bucket <- Sys.getenv("WORKSPACE_BUCKET")

for (fname in c("dt_agg_pums_v2.csv", "dt_agg_aou_v3.csv", "aou_raking_dt.csv")) {
  system(paste0("gsutil cp ", my_bucket, "/data/", fname, " ."), intern = TRUE)
}

dt_agg_pums <- read_csv("dt_agg_pums_v2.csv")
dt_agg_aou  <- read_csv("dt_agg_aou_v3.csv")
dt_raw_aou  <- read_csv("aou_raking_dt.csv")

########################################################################
## Apply factor levels (same harmonization as the Global RAILS run)
########################################################################

harmonize_factors <- function(df) {
  df %>%
    mutate(
      sex      = factor(sex,      levels = c("Female", "Male")),
      race_eth = factor(race_eth, levels = c("Hispanic", "NH Asian", "NH Black", "NH White", "Others")),
      income   = factor(income,   levels = c("<35k", "35k-50k", "50k-75k", "75k-100k", ">100k")),
      agegroup = factor(agegroup, levels = c("18-24", "25-44", "45-64", "65-74", "75+")),
      edu      = factor(edu,      levels = c("Less than highschool", "Some highschool",
                                             "Highschool graduate", "Some college",
                                             "College graduate or advanced")),
      homeown  = factor(ifelse(homeown == "Other", "Others", homeown),
                        levels = c("Own", "Rent", "Others")),
      region   = factor(region,   levels = c("Northeast", "Midwest", "South", "West"))
    ) %>%
    na.omit()
}

dt_agg_pums <- harmonize_factors(dt_agg_pums)
dt_agg_aou  <- harmonize_factors(dt_agg_aou)

########################################################################
## Individual-level AoU (original cleaning chain, for the final join)
########################################################################

south     <- c("AL","AR","FL","GA","KY","LA","MS","NC","SC","TN","TX","VA","WV")
midwest   <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI")
northeast <- c("CT","DE","ME","MD","MA","NH","NJ","NY","PA","RI","VT")
west      <- c("AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY")
state_region_map <- c(
  setNames(rep("South",     length(south)),     south),
  setNames(rep("Midwest",   length(midwest)),   midwest),
  setNames(rep("Northeast", length(northeast)), northeast),
  setNames(rep("West",      length(west)),      west)
)

dt_aou <- dt_raw_aou %>%
  select(!c(gender, careplace)) %>%
  na.omit() %>%
  mutate(
    race_eth = case_when(
      ethnicity == "Yes"                       ~ "Hispanic",
      race      == "Asian"                     ~ "NH Asian",
      race      == "Black or African American" ~ "NH Black",
      race      == "White"                     ~ "NH White",
      TRUE                                     ~ "Others"
    ),
    agegroup = case_when(
      age <= 24            ~ "18-24",
      age > 24 & age <= 44 ~ "25-44",
      age > 44 & age <= 64 ~ "45-64",
      age > 64 & age <= 74 ~ "65-74",
      TRUE                 ~ "75+"
    ),
    region = state_region_map[state]
  ) %>%
  harmonize_factors()

########################################################################
## Restrict to the chosen subgroup
## SUB_VAR drops out of the model (constant within the stratum); all
## remaining covariates keep their full factor levels.
########################################################################

names_univar <- setdiff(names_univar_all, SUB_VAR)

dt_agg_pums_sub <- dt_agg_pums %>% filter(.data[[SUB_VAR]] == SUB_LEVEL)
dt_agg_aou_sub  <- dt_agg_aou  %>% filter(.data[[SUB_VAR]] == SUB_LEVEL)
dt_aou_sub      <- dt_aou      %>% filter(.data[[SUB_VAR]] == SUB_LEVEL)

stopifnot(nrow(dt_agg_pums_sub) > 0, nrow(dt_agg_aou_sub) > 0)

## Subgroup population total — the nsiz of this analysis
nsiz_sub <- sum(dt_agg_pums_sub$weight)
message(SUB_VAR, " = ", SUB_LEVEL, ": ", nrow(dt_agg_aou_sub), " AoU cells, ",
        sum(dt_agg_aou_sub$weight), " participants, nsiz = ", format(nsiz_sub, big.mark = ","))

########################################################################
## Precompute all subgroup PUMS margins (one-way + two-way + three-way)
########################################################################

twovars   <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
threevars <- combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))

max_formula <- formula(paste0("~", paste(c(names_univar, twovars, threevars), collapse = "+")))
mat_max     <- sparse.model.matrix(max_formula, data = dt_agg_pums_sub, keep.order = TRUE)
pop_totals_sub <- setNames(
  as.numeric(Matrix::crossprod(mat_max, dt_agg_pums_sub$weight)),
  colnames(mat_max)
)

########################################################################
## Run subgroup RAILS
########################################################################

result_sub <- fun.rails.threeway(
  dt_agg_aou   = dt_agg_aou_sub,
  dt_agg_pums  = dt_agg_pums_sub,
  pop_totals   = pop_totals_sub,
  names_univar = names_univar,
  alpha        = 0.05,
  nsiz         = nsiz_sub
)

## Diagnostics: total ≈ nsiz_sub, and the model actually calibrated
sum(result_sub$d_rails * result_sub$weight, na.rm = TRUE)
result_sub$calibrated_terms[1]
fun.out(result_sub$d_rails)

########################################################################
## Attach weights to individual-level subgroup AoU data
## All d_* columns are already per-individual.
########################################################################

dt_out_sub <- dt_aou_sub %>%
  left_join(
    result_sub %>%
      rename(
        w_unweighted = d_unweighted,
        w_cal1       = d_cal1,
        w_cal2       = d_cal2,
        w_nps1       = d_nps1,
        w_nps2       = d_nps2,
        w_nps1_rake  = d_nps1_rake,
        w_nps2_rake  = d_nps2_rake,
        w_rails      = d_rails
      ) %>%
      select(all_of(names_univar), w_unweighted, w_cal1, w_cal2,
             w_nps1, w_nps2, w_nps1_rake, w_nps2_rake, w_rails,
             selected_terms, calibrated_terms),
    by = names_univar
  )

########################################################################
## Save results to bucket
########################################################################

write_excel_csv(dt_out_sub, OUT_FILE)
system(paste0("gsutil cp ./", OUT_FILE, " ", my_bucket, "/data/"), intern = TRUE)
system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = TRUE)
