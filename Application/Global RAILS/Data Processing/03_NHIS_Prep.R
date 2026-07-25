#### 03_NHIS_Prep.R — Download and prepare 2020 NHIS reference data (OPTIONAL)
#### Needed only for the hybrid design: NHIS as the propensity reference
#### sample while raking calibrates to PUMS population totals.
#### Can be run locally or inside the AoU Researcher Workbench.

library(tidyverse)
library(survey)

########################################################################
## Download 2020 NHIS adult file from NCHS (cached after first run)
########################################################################

if (!"NHIS_2020_raw.csv" %in% list.files(getwd(), recursive = TRUE)) {

  zip_url   <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NHIS/2020/adult20csv.zip"
  temp      <- tempfile()
  download.file(zip_url, temp)
  unzip_dir <- tempdir()
  unzip(temp, exdir = unzip_dir)
  csv_file  <- list.files(unzip_dir, pattern = "20.csv$", full.names = TRUE)
  df        <- read.csv(csv_file)
  write.csv(df, "NHIS_2020_raw.csv", row.names = FALSE)

} else {
  df <- read.csv("NHIS_2020_raw.csv")
}

########################################################################
## Select, rename, and recode variables to the shared categories
########################################################################

dt <- df %>%
  select(WTFA_A, AGEP_A, SEX_A, EDUC_A, HISP_A, HISPALLP_A, PHSTAT_A,
         REGION, INCGRP_A, USPLKIND_A, HOUTENURE_A) %>%
  rename(
    weight        = WTFA_A,
    age           = AGEP_A,
    sex           = SEX_A,
    edu           = EDUC_A,
    hispanic_ind  = HISP_A,
    race_eth      = HISPALLP_A,
    health_status = PHSTAT_A,
    region        = REGION,
    income        = INCGRP_A,
    careplace     = USPLKIND_A,
    homeown       = HOUTENURE_A
  ) %>%
  mutate(sex = case_when(
    sex == 1 ~ "Male",
    sex == 2 ~ "Female",
    TRUE     ~ NA_character_
  )) %>%
  mutate(edu = case_when(
    edu %in% 8:11 ~ "College graduate or advanced",
    edu %in% 5:7  ~ "Some college",
    edu %in% 3:4  ~ "Highschool graduate",
    edu == 2      ~ "Some highschool",
    edu %in% 0:1  ~ "Less than highschool",
    TRUE          ~ NA_character_
  )) %>%
  mutate(race_eth = case_when(
    race_eth == 1      ~ "Hispanic",
    race_eth == 2      ~ "NH White",
    race_eth == 3      ~ "NH Black",
    race_eth == 4      ~ "NH Asian",
    race_eth %in% 5:7  ~ "Others",
    TRUE               ~ NA_character_
  )) %>%
  mutate(health_status = case_when(
    health_status == 1 ~ "Excellent",
    health_status == 2 ~ "Very Good",
    health_status == 3 ~ "Good",
    health_status == 4 ~ "Fair",
    health_status == 5 ~ "Poor",
    TRUE               ~ NA_character_
  )) %>%
  mutate(region = case_when(
    region == 1 ~ "Northeast",
    region == 2 ~ "Midwest",
    region == 3 ~ "South",
    region == 4 ~ "West",
    TRUE        ~ NA_character_
  )) %>%
  mutate(income = case_when(
    income == 1 ~ "<35k",
    income == 2 ~ "35k-50k",
    income == 3 ~ "50k-75k",
    income == 4 ~ "75k-100k",
    income == 5 ~ ">100k",
    TRUE        ~ NA_character_
  )) %>%
  mutate(careplace = case_when(
    careplace == 1     ~ "Doctors office",
    careplace == 2     ~ "Urgent",
    careplace == 3     ~ "Emergency",
    careplace == 6     ~ "None",
    careplace %in% 4:5 ~ "Others",
    TRUE               ~ NA_character_
  )) %>%
  mutate(homeown = case_when(
    homeown == 1 ~ "Own",
    homeown == 2 ~ "Rent",
    homeown == 3 ~ "Others",
    TRUE         ~ NA_character_
  ))

########################################################################
## Save recoded individual-level NHIS data
########################################################################

write_excel_csv(dt, "NHIS_2020.csv")

########################################################################
## Analysis dataset: drop unused variables, age groups, factor levels
########################################################################

dt_nhis <- dt %>%
  select(!c(hispanic_ind, health_status, careplace)) %>%
  na.omit() %>%
  mutate(agegroup = case_when(
    age <= 24            ~ "18-24",
    age > 24 & age <= 44 ~ "25-44",
    age > 44 & age <= 64 ~ "45-64",
    age > 64 & age <= 74 ~ "65-74",
    TRUE                 ~ "75+"
  )) %>%
  mutate(agegroup = factor(agegroup, levels = c("18-24","25-44","45-64","65-74","75+"))) %>%
  mutate(sex      = factor(sex,      levels = c("Female","Male"))) %>%
  mutate(race_eth = factor(race_eth, levels = c("Hispanic","NH Asian","NH Black","NH White","Others"))) %>%
  mutate(income   = factor(income,   levels = c("<35k","35k-50k","50k-75k","75k-100k",">100k"))) %>%
  mutate(edu      = factor(edu,      levels = c("Less than highschool","Some highschool",
                                                "Highschool graduate","Some college",
                                                "College graduate or advanced"))) %>%
  mutate(homeown  = factor(homeown,  levels = c("Own","Rent","Others"))) %>%
  mutate(region   = factor(region,   levels = c("Northeast","Midwest","South","West"))) %>%
  na.omit()

########################################################################
## Aggregate to weighted cell counts (the dt_s role in fun.rails.threeway)
########################################################################

names_univar <- c("agegroup","sex","edu","homeown","income","race_eth","region")
cat_formula  <- formula(paste0("weight ~", paste(names_univar, collapse = "+")))
dt_agg_nhis  <- aggregate(cat_formula, data = dt_nhis, sum)

write_excel_csv(dt_agg_nhis, "dt_agg_nhis_v2.csv")

########################################################################
## Save to Workbench bucket (skipped when running locally)
########################################################################

my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
if (nchar(my_bucket) > 0) {
  system(paste0("gsutil cp ./NHIS_2020.csv ",      my_bucket, "/data/"), intern = TRUE)
  system(paste0("gsutil cp ./dt_agg_nhis_v2.csv ", my_bucket, "/data/"), intern = TRUE)
  system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = TRUE)
} else {
  message("WORKSPACE_BUCKET not set — skipping upload. Copy NHIS_2020.csv and dt_agg_nhis_v2.csv to the Workbench manually.")
}
