#### 04_Global_RAILS.R — Global RAILS Analysis
#### Run inside the AoU Researcher Workbench

library(tidyverse)
library(survey)
library(Matrix)
library(plotly)
source("AoU_Fun.R")

########################################################################
## Load aggregated datasets from Workspace bucket
########################################################################

my_bucket <- Sys.getenv("WORKSPACE_BUCKET")

for (fname in c("dt_agg_pums_v2.csv", "dt_agg_aou_v3.csv", "aou_raking_dt.csv")) {
  system(paste0("gsutil cp ", my_bucket, "/data/", fname, " ."), intern = TRUE)
}

dt_agg_pums <- read_csv("dt_agg_pums_v2.csv")
dt_agg_aou  <- read_csv("dt_agg_aou_v3.csv")
dt_raw_aou  <- read_csv("aou_raking_dt.csv")

########################################################################
## Apply factor levels
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

names_univar <- c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region")

########################################################################
## Prepare individual-level AoU data
########################################################################

## Define the states in each region
south     <- c("AL","AR","FL","GA","KY","LA","MS","NC","SC","TN","TX","VA","WV")
midwest   <- c("IL","IN","IA","KS","MI","MN","MO","NE","ND","OH","SD","WI")
northeast <- c("CT","DE","ME","MD","MA","NH","NJ","NY","PA","RI","VT")
west      <- c("AK","AZ","CA","CO","HI","ID","MT","NV","NM","OR","UT","WA","WY")

## Create a named vector to map states to regions
state_region_map <- c(
  setNames(rep("South",     length(south)),     south),
  setNames(rep("Midwest",   length(midwest)),   midwest),
  setNames(rep("Northeast", length(northeast)), northeast),
  setNames(rep("West",      length(west)),      west)
)

dt_aou <- dt_raw_aou %>%
  select(!c(gender, careplace)) %>%
  na.omit()

dt_aou <- dt_aou %>%
  mutate(race_eth = case_when(
    ethnicity == "Yes"                       ~ "Hispanic",
    race      == "Asian"                     ~ "NH Asian",
    race      == "Black or African American" ~ "NH Black",
    race      == "White"                     ~ "NH White",
    TRUE                                     ~ "Others"
  ))

dt_aou <- dt_aou %>%
  mutate(agegroup = case_when(
    age <= 24            ~ "18-24",
    age > 24 & age <= 44 ~ "25-44",
    age > 44 & age <= 64 ~ "45-64",
    age > 64 & age <= 74 ~ "65-74",
    TRUE                 ~ "75+"
  ))

dt_aou <- dt_aou %>%
  mutate(sex      = factor(sex,      levels = c("Female","Male"))) %>%
  mutate(race_eth = factor(race_eth, levels = c("Hispanic","NH Asian","NH Black","NH White","Others"))) %>%
  mutate(income   = factor(income,   levels = c("<35k","35k-50k","50k-75k","75k-100k",">100k"))) %>%
  mutate(agegroup = factor(agegroup, levels = c("18-24","25-44","45-64","65-74","75+"))) %>%
  mutate(edu      = factor(edu,      levels = c("Less than highschool","Some highschool",
                                                "Highschool graduate","Some college",
                                                "College graduate or advanced"))) %>%
  mutate(homeown  = ifelse(homeown == "Other", "Others", homeown)) %>%
  mutate(homeown  = factor(homeown,  levels = c("Own","Rent","Others"))) %>%
  na.omit()

dt_aou <- dt_aou %>%
  mutate(region = state_region_map[state]) %>%
  mutate(region = factor(region, levels = c("Northeast","Midwest","South","West"))) %>%
  na.omit() %>%
  mutate(weight = 1)

########################################################################
## Precompute all population margins from PUMS (one-way + two-way + three-way)
## This is done once here so fun.rails.threeway never needs sparse.model.matrix internally.
########################################################################

twovars   <- combn(names_univar, 2, FUN = function(x) paste(x, collapse = ":"))
threevars <- combn(names_univar, 3, FUN = function(x) paste(x, collapse = ":"))

max_formula <- formula(paste0("~", paste(c(names_univar, twovars, threevars), collapse = "+")))
mat_max     <- sparse.model.matrix(max_formula, data = dt_agg_pums, keep.order = TRUE)
pop_totals  <- as.numeric(Matrix::crossprod(mat_max, dt_agg_pums$weight))
names(pop_totals) <- colnames(mat_max)

########################################################################
## Run Global RAILS
########################################################################

result_rails <- fun.rails.threeway(
  dt_agg_aou   = dt_agg_aou,
  dt_agg_pums  = dt_agg_pums,
  pop_totals   = pop_totals,
  names_univar = names_univar,
  alpha        = 0.05
)

########################################################################
## Attach RAILS + benchmark weights to individual-level AoU data
########################################################################

## result_rails has one row per covariate cell; all d_* columns are already
## PER-INDIVIDUAL weights (fun.rails.threeway divides cell totals by cell counts).
## Benchmarks: cal-1/cal-2 (raking only), nps-1/nps-2 (PS only),
## nps-1-rake/nps-2-rake (PS + raking), and w_rails (RAILS).
dt_aou_weighted <- dt_aou %>%
  left_join(
    result_rails %>%
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

write_excel_csv(dt_aou_weighted, "global_rails_weights.csv")
system(paste0("gsutil cp ./global_rails_weights.csv ", my_bucket, "/data/"), intern = TRUE)
system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = TRUE)

########################################################################
## OPTIONAL — Hybrid design: NHIS as propensity reference, PUMS totals
## Requires dt_agg_nhis_v2.csv in the bucket (see 03_NHIS_Prep.R).
## Set run_hybrid <- TRUE to enable.
########################################################################

run_hybrid <- FALSE

if (run_hybrid) {

  system(paste0("gsutil cp ", my_bucket, "/data/dt_agg_nhis_v2.csv ."), intern = TRUE)
  dt_agg_nhis <- read_csv("dt_agg_nhis_v2.csv") %>% harmonize_factors()

  ## Population total: same definition as the PUMS-only run — the NPS weights
  ## are scaled to the PUMS adult total so the scale matches pop_totals,
  ## rather than the NHIS weight sum.
  nsiz_hybrid <- sum(dt_agg_pums$weight)

  result_hybrid <- fun.rails.threeway(
    dt_agg_aou   = dt_agg_aou,
    dt_agg_pums  = dt_agg_nhis,   # NHIS cells in the reference-sample slot
    pop_totals   = pop_totals,    # raking targets computed from PUMS
    names_univar = names_univar,
    alpha        = 0.05,
    nsiz         = nsiz_hybrid
  )

  ## d_* columns are already per-individual weights
  dt_aou_weighted_hybrid <- dt_aou %>%
    left_join(
      result_hybrid %>%
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

  write_excel_csv(dt_aou_weighted_hybrid, "hybrid_rails_weights.csv")
  system(paste0("gsutil cp ./hybrid_rails_weights.csv ", my_bucket, "/data/"), intern = TRUE)
}

