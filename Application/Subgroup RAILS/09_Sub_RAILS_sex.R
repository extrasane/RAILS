#### 09_Sub_RAILS_sex.R — Female-only RAILS analysis (breast cancer study)
#### Run inside the AoU Researcher Workbench.
#### Subgroup variable: sex. Restricted to the FEMALE stratum, so the wrapper
#### fun.sub.rails.threeway runs exactly one level.
####
#### Requires Sub_AoU_Fun.R AND AoU_Fun.R in this directory
#### (Sub_AoU_Fun.R sources AoU_Fun.R).

library(tidyverse)
library(survey)
library(Matrix)
source("Sub_AoU_Fun.R")          # fun.sub.rails.threeway (+ AoU_Fun.R)

########################################################################
## Load aggregated datasets from Workspace bucket
########################################################################

my_bucket <- Sys.getenv("WORKSPACE_BUCKET")

for (fname in c("dt_agg_pums_v2.csv", "dt_agg_aou_v2.csv", "aou_raking_dt.csv")) {
  system(paste0("gsutil cp ", my_bucket, "/data/", fname, " ."), intern = TRUE)
}

dt_agg_pums <- read_csv("dt_agg_pums_v2.csv")
dt_agg_aou  <- read_csv("dt_agg_aou_v2.csv")
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
## Restrict to the FEMALE population
## sex is the subgroup variable (single level here); the wrapper excludes it
## from the model automatically.
########################################################################

names_univar <- c("agegroup", "edu", "homeown", "income", "race_eth", "region")

dt_agg_pums_fem <- dt_agg_pums %>% filter(sex == "Female")
dt_agg_aou_fem  <- dt_agg_aou  %>% filter(sex == "Female")
dt_aou_fem      <- dt_aou      %>% filter(sex == "Female")

########################################################################
## Run Female RAILS via the subgroup wrapper (subgroup_var = "sex")
## Pre-filtered to Female, so the wrapper loops one level and builds the
## female PUMS margins + nsiz internally.
########################################################################

result_sex <- fun.sub.rails.threeway(
  dt_agg_aou   = dt_agg_aou_fem,
  dt_agg_pums  = dt_agg_pums_fem,
  subgroup_var = "sex",
  names_univar = names_univar,
  alpha        = 0.05
)

## Diagnostics: per-subgroup totals and the calibrated model
sum(result_sex$d_rails * result_sex$weight, na.rm = TRUE)
result_sex$calibrated_terms[1]
fun.out(result_sex$d_rails)

########################################################################
## Attach weights to individual-level FEMALE AoU data
## All d_* columns are already per-individual.
########################################################################

dt_out_sex <- dt_aou_fem %>%
  left_join(
    result_sex %>%
      rename(
        w_unweighted = d_unweighted,
        w_cal1       = d_cal1,
        w_cal2       = d_cal2,
        w_nps1       = d_nps1,
        w_nps2       = d_nps2,
        w_nps1_rake  = d_nps1_rake,
        w_nps2_rake  = d_nps2_rake,
        w_subrails   = d_rails
      ) %>%
      select(all_of(names_univar), w_unweighted, w_cal1, w_cal2,
             w_nps1, w_nps2, w_nps1_rake, w_nps2_rake, w_subrails,
             selected_terms, calibrated_terms),
    by = names_univar
  )

########################################################################
## Save results — LOCAL only for now.
## Bucket upload annotated out until the input data version is confirmed.
########################################################################

write_excel_csv(dt_out_sex, "../data/dt_sub_aou_sex_female.csv")

# system(paste0("gsutil cp ../data/dt_sub_aou_sex_female.csv ", my_bucket, "/data/"), intern = TRUE)
# system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = TRUE)
