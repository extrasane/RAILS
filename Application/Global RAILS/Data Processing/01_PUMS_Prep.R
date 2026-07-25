#### 01_PUMS_Prep.R — Download and prepare 2022 ACS PUMS reference data
#### Can be run locally or in any R environment (does not require the Workbench)

# install.packages(c("httr", "jsonlite", "dplyr", "survey"))
# NOTE: if a cached PUMS_2022.csv exists from an older extraction WITHOUT the
# REGION column, delete it so the download re-runs with REGION included.
library(httr)
library(jsonlite)
library(dplyr)
library(tidyverse)
library(Matrix)
library(survey)

########################################################################
## Download PUMS from Census API (cached after first run)
########################################################################

# Request your own free key at: https://api.census.gov/data/key_signup.html
CENSUS_KEY <- "YOUR_CENSUS_API_KEY"

if (!"PUMS_2022.csv" %in% list.files(getwd(), recursive = TRUE)) {

  API_URL <- paste0(
    "https://api.census.gov/data/2022/acs/acs1/pums?",
    "get=PWGTP,HINCP,AGEP,RACSOR,RACAIAN,RACASN,RACBLK,RACWHT,",
    "TEN,SEX,RACNH,HISP,SCHL,RACPI,REGION",
    "&ucgid=0200000US1,0200000US2,0200000US3,0200000US4",
    "&key=", CENSUS_KEY
  )

  raw_data <- GET(API_URL)
  stop_for_status(raw_data)

  dt_list <- fromJSON(rawToChar(raw_data$content))
  dt_list <- data.frame(dt_list)
  temp    <- dt_list[1, ]
  dt_list <- dt_list[-1, ]
  colnames(dt_list) <- temp

  write.csv(dt_list, "PUMS_2022.csv", row.names = FALSE)

} else {
  dt_list <- read.csv("PUMS_2022.csv")
}

## The API response appends geography columns (e.g. lowercase "region" or
## "ucgid") alongside the requested REGION variable. Keep only the first of
## any case-insensitive duplicate so rename(region = REGION) has one target.
dt_list <- dt_list[, !duplicated(toupper(colnames(dt_list)))]

########################################################################
## Recode and harmonize variables
##
## DISCLAIMER: All variable codings below (SCHL, TEN, SEX, HISP, RAC*,
## HINCP, AGEP, ST) follow the 2022 ACS 1-year PUMS Data Dictionary:
##   https://www.census.gov/programs-surveys/acs/microdata/documentation.html
## If you use a different PUMS year, verify each code against that
## year's dictionary — category codes can change between releases.
########################################################################

dt <- dt_list %>%
  rename(
    w_pums               = PWGTP,
    income               = HINCP,
    age                  = AGEP,
    race_other           = RACSOR,
    race_American_Indian = RACAIAN,
    race_Asian           = RACASN,
    race_Black           = RACBLK,
    race_White           = RACWHT,
    sex                  = SEX,
    race_hawaiian        = RACNH,
    hispanic             = HISP,
    edu                  = SCHL,
    race_Pacific_Islander = RACPI,
    region               = REGION,
    homeown              = TEN
  ) %>%
  ## Coerce coded fields to numeric BEFORE any comparison. The API returns
  ## everything as character with zero-padded codes (e.g. HISP "01"), so a
  ## character comparison like hispanic != 1 silently misclassifies — the
  ## cached read.csv path auto-converts to integer, but the fresh-download
  ## path does not. Explicit coercion makes both paths identical.
  mutate(across(c(age, hispanic, sex, race_White, race_Black, race_Asian,
                  income, edu, region, homeown, w_pums), as.numeric)) %>%
  filter(age > 17) %>%
  mutate(eth = case_when(
    hispanic != 1 ~ "Hispanic",
    TRUE          ~ "Non-Hispanic"
  )) %>%
  mutate(race = case_when(
    race_White == 1 ~ "White",
    race_Black == 1 ~ "Black",
    race_Asian == 1 ~ "Asian",
    TRUE            ~ "Others"
  )) %>%
  mutate(sex = case_when(
    sex == 1 ~ "Male",
    sex == 2 ~ "Female",
    TRUE     ~ NA_character_
  )) %>%
  mutate(race_eth = case_when(
    eth == "Hispanic"                          ~ "Hispanic",
    race == "White" & eth == "Non-Hispanic"    ~ "NH White",
    race == "Black" & eth == "Non-Hispanic"    ~ "NH Black",
    race == "Asian" & eth == "Non-Hispanic"    ~ "NH Asian",
    TRUE                                       ~ "Others"
  )) %>%
  mutate(income = as.numeric(income)) %>%
  mutate(income = case_when(
    income >= -60000 & income < 35000  ~ "<35k",
    income >= 35000  & income < 50000  ~ "35k-50k",
    income >= 50000  & income < 75000  ~ "50k-75k",
    income >= 75000  & income < 100000 ~ "75k-100k",
    income >= 100000                   ~ ">100k",
    TRUE                               ~ NA_character_
  )) %>%
  mutate(edu = as.numeric(edu)) %>%
  mutate(edu = case_when(
    edu < 12                    ~ "Less than highschool",
    edu >= 12 & edu < 16        ~ "Some highschool",
    edu == 16 | edu == 17       ~ "Highschool graduate",
    edu >= 18 & edu < 21        ~ "Some college",
    edu >= 21                   ~ "College graduate or advanced",
    TRUE                        ~ NA_character_
  )) %>%
  mutate(homeown = as.numeric(homeown)) %>%
  mutate(homeown = case_when(
    homeown == 1 | homeown == 2 ~ "Own",
    homeown == 3                ~ "Rent",
    homeown == 4                ~ "Others",
    TRUE                        ~ NA_character_
  )) %>%
  ## Official Census region from the PUMS REGION variable
  mutate(region = as.numeric(region)) %>%
  mutate(region = case_when(
    region == 1 ~ "Northeast",
    region == 2 ~ "Midwest",
    region == 3 ~ "South",
    region == 4 ~ "West",
    TRUE        ~ NA_character_
  )) %>%
  mutate(age = as.numeric(age))

## Adult population total, computed BEFORE na.omit — the anchor the weights
## are rescaled to, and the nsiz used downstream (matches the original
## extraction, which produced 261,048,645 for 2022)
n_adult <- sum(as.numeric(dt$w_pums)[dt$age >= 18], na.rm = TRUE)

########################################################################
## Create age groups, apply factor levels, drop incomplete rows
## (dt is adults-only: age > 17 filtered at the top of the recode chain)
########################################################################

dt_pums <- dt %>%
  mutate(agegroup = case_when(
    age >= 18 & age <= 24  ~ "18-24",
    age > 24 & age <= 44   ~ "25-44",
    age > 44 & age <= 64   ~ "45-64",
    age > 64 & age <= 74   ~ "65-74",
    TRUE                   ~ "75+"
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
  mutate(w_pums   = as.numeric(w_pums)) %>%
  na.omit()

## Rescale ADULT weights to the pre-na.omit adult total (original anchor)
dt_pums <- dt_pums %>%
  mutate(w_pums = w_pums / sum(w_pums) * n_adult)

########################################################################
## Aggregate to weighted cell counts and save to Workbench bucket
########################################################################

dt_pums_analysis <- dt_pums %>%
  rename(weight = w_pums) %>%
  select(weight, agegroup, sex, edu, homeown, income, race_eth, region)

names_univar <- c("agegroup","sex","edu","homeown","income","race_eth","region")
cat_formula  <- formula(paste0("weight ~", paste(names_univar, collapse = "+")))
dt_agg_pums  <- aggregate(cat_formula, data = dt_pums_analysis, sum)

write_excel_csv(dt_agg_pums, "dt_agg_pums_v2.csv")

my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
if (nchar(my_bucket) > 0) {
  system(paste0("gsutil cp ./dt_agg_pums_v2.csv ", my_bucket, "/data/"), intern = TRUE)
  system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = TRUE)
} else {
  message("WORKSPACE_BUCKET not set — skipping upload. Copy dt_agg_pums_v2.csv to the Workbench manually.")
}
