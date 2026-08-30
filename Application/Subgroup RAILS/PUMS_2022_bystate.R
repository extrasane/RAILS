####PUMS_Prep_bystate.R — 2022 ACS PUMS reference data, BY STATE
#### Can be run locally or in any R environment (does not require the Workbench).
####
#### Purpose: produce ../data/PUMS_2022_bystate.csv, the per-person state
#### population file consumed by 10_State_Prevalence_Maps.R to build the
#### state population shares rho_a = N^a / N (STATE_POP_WEIGHT = "PWGTP").
####
#### This mirrors the recoding grammar of
####   Application/Global RAILS/Data Processing/01_PUMS_Prep.R
#### with three deliberate differences for the state-level use case:
####   1. Geography is pulled at STATE level (ucgid 0400000US..), and the
####      state variable ST is requested and mapped to a 2-letter USPS code
####      that matches the AoU `state` column the maps join on.
####   2. The person-weight column is kept as PWGTP (NOT renamed to w_pums),
####      so 10_State_Prevalence_Maps.R finds it via STATE_POP_WEIGHT.
####   3. No na.omit() before writing: rho_a is a *population* share, so the
####      raw ACS person weights must cover ALL adults in each state. Dropping
####      rows with a missing demographic (income/edu/homeown) would undercount
####      state populations. Demographic columns are still recoded and kept for
####      cross-checks, and may contain NA.

# install.packages(c("httr", "jsonlite", "dplyr", "tidyverse"))
library(httr)
library(jsonlite)
library(dplyr)
library(tidyverse)

########################################################################
## State FIPS -> USPS lookup (50 states + DC), self-contained so this
## script needs no tigris/network state table. Also drives the ucgid list.
########################################################################

fips2usps <- c(
  "01"="AL","02"="AK","04"="AZ","05"="AR","06"="CA","08"="CO","09"="CT",
  "10"="DE","11"="DC","12"="FL","13"="GA","15"="HI","16"="ID","17"="IL",
  "18"="IN","19"="IA","20"="KS","21"="KY","22"="LA","23"="ME","24"="MD",
  "25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT","31"="NE",
  "32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND",
  "39"="OH","40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD",
  "47"="TN","48"="TX","49"="UT","50"="VT","51"="VA","53"="WA","54"="WV",
  "55"="WI","56"="WY")

state_ucgids <- paste0("0400000US", names(fips2usps))

########################################################################
## Download PUMS from Census API (cached after first run)
########################################################################

# Request your own free key at: https://api.census.gov/data/key_signup.html
# CENSUS_KEY <- YOURKEY

PUMS_CACHE <- "../data/PUMS_2022_state.csv"

## Fetch one state's PUMS records. We query ONE geography per call: the ACS
## API caps how many `ucgid` values a single request accepts, so asking for all
## 51 states at once makes the gateway return an HTML error page (which then
## fails to parse as JSON). Requesting each state separately stays under the cap
## and lets us surface the real server message if a call fails.
GET_VARS <- paste0(
  "PWGTP,HINCP,AGEP,RACSOR,RACAIAN,RACASN,RACBLK,RACWHT,",
  "TEN,SEX,RACNH,HISP,SCHL,RACPI,REGION,ST"              # ST = state of residence
)

fetch_pums_state <- function(ucgid) {
  key_q <- if (nzchar(CENSUS_KEY) && CENSUS_KEY != "YOUR_CENSUS_API_KEY")
    paste0("&key=", CENSUS_KEY) else ""
  url  <- paste0("https://api.census.gov/data/2022/acs/acs1/pums?get=",
                 GET_VARS, "&ucgid=", ucgid, key_q)
  resp <- GET(url)
  txt  <- content(resp, as = "text", encoding = "UTF-8")
  ## A valid PUMS response is a JSON array, so it must start with '['. Anything
  ## else (HTML, a plain-text complaint) is an error we should show verbatim.
  if (http_error(resp) || !startsWith(trimws(txt), "[")) {
    stop("Census API error for ucgid=", ucgid, " (HTTP ", status_code(resp), ").\n",
         "First 300 chars of response:\n", substr(trimws(txt), 1, 300),
         call. = FALSE)
  }
  m  <- fromJSON(txt)                                    # character matrix
  df <- as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE)
  colnames(df) <- m[1, ]
  ## A state-level ucgid query returns ST twice: once as the requested variable
  ## and once as the appended geography column (same value). Drop the duplicate
  ## HERE, before bind_rows repairs the names to ST...16 / ST...17.
  df[, !duplicated(toupper(colnames(df))), drop = FALSE]
}

if (!file.exists(PUMS_CACHE)) {
  message("Downloading PUMS 2022 by state (", length(state_ucgids), " calls)...")
  parts <- lapply(seq_along(state_ucgids), function(i) {
    message("  [", i, "/", length(state_ucgids), "] ",
            names(fips2usps)[i], " (", fips2usps[i], ")")
    fetch_pums_state(state_ucgids[i])
  })
  dt_list <- dplyr::bind_rows(parts)
  write.csv(dt_list, PUMS_CACHE, row.names = FALSE)
} else {
  dt_list <- read.csv(PUMS_CACHE)
}

## Safety net (also repairs a stale cache written before the fetch-side de-dup):
## strip any "...N" name-repair suffix that bind_rows/read.csv may have added to
## duplicate geography columns (e.g. ST...16 / ST...17), then keep the first of
## each case-insensitive duplicate so rename() has a single target.
colnames(dt_list) <- sub("\\.\\.\\..*$", "", colnames(dt_list))
dt_list <- dt_list[, !duplicated(toupper(colnames(dt_list)))]

########################################################################
## Recode and harmonize variables (2022 ACS 1-year PUMS Data Dictionary)
## Codings follow Global RAILS 01_PUMS_Prep.R exactly.
########################################################################

dt <- dt_list %>%
  rename(
    income                = HINCP,
    age                   = AGEP,
    race_other            = RACSOR,
    race_American_Indian  = RACAIAN,
    race_Asian            = RACASN,
    race_Black            = RACBLK,
    race_White            = RACWHT,
    sex                   = SEX,
    race_hawaiian         = RACNH,
    hispanic              = HISP,
    edu                   = SCHL,
    race_Pacific_Islander = RACPI,
    region                = REGION,
    homeown               = TEN,
    st_fips               = ST
    ## NOTE: PWGTP is intentionally NOT renamed — the maps script reads it by name.
  ) %>%
  ## Coerce coded fields to numeric BEFORE any comparison (fresh-download path
  ## returns zero-padded character codes; cached read.csv path returns integers).
  mutate(across(c(age, hispanic, sex, race_White, race_Black, race_Asian,
                  income, edu, region, homeown, st_fips, PWGTP), as.numeric)) %>%
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
    eth == "Hispanic"                       ~ "Hispanic",
    race == "White" & eth == "Non-Hispanic" ~ "NH White",
    race == "Black" & eth == "Non-Hispanic" ~ "NH Black",
    race == "Asian" & eth == "Non-Hispanic" ~ "NH Asian",
    TRUE                                    ~ "Others"
  )) %>%
  mutate(income = case_when(
    income >= -60000 & income < 35000  ~ "<35k",
    income >= 35000  & income < 50000  ~ "35k-50k",
    income >= 50000  & income < 75000  ~ "50k-75k",
    income >= 75000  & income < 100000 ~ "75k-100k",
    income >= 100000                   ~ ">100k",
    TRUE                               ~ NA_character_
  )) %>%
  mutate(edu = case_when(
    edu < 12              ~ "Less than highschool",
    edu >= 12 & edu < 16  ~ "Some highschool",
    edu == 16 | edu == 17 ~ "Highschool graduate",
    edu >= 18 & edu < 21  ~ "Some college",
    edu >= 21             ~ "College graduate or advanced",
    TRUE                  ~ NA_character_
  )) %>%
  mutate(homeown = case_when(
    homeown == 1 | homeown == 2 ~ "Own",
    homeown == 3                ~ "Rent",
    homeown == 4                ~ "Others",
    TRUE                        ~ NA_character_
  )) %>%
  ## Official Census region from the PUMS REGION variable (not a hand-built map)
  mutate(region = case_when(
    region == 1 ~ "Northeast",
    region == 2 ~ "Midwest",
    region == 3 ~ "South",
    region == 4 ~ "West",
    TRUE        ~ NA_character_
  )) %>%
  ## State of residence: FIPS -> 2-letter USPS (matches the AoU `state` column)
  mutate(state = unname(fips2usps[sprintf("%02d", st_fips)])) %>%
  mutate(agegroup = case_when(
    age >= 18 & age <= 24 ~ "18-24",
    age > 24 & age <= 44  ~ "25-44",
    age > 44 & age <= 64  ~ "45-64",
    age > 64 & age <= 74  ~ "65-74",
    TRUE                  ~ "75+"
  )) %>%
  mutate(
    agegroup = factor(agegroup, levels = c("18-24","25-44","45-64","65-74","75+")),
    sex      = factor(sex,      levels = c("Female","Male")),
    race_eth = factor(race_eth, levels = c("Hispanic","NH Asian","NH Black","NH White","Others")),
    income   = factor(income,   levels = c("<35k","35k-50k","50k-75k","75k-100k",">100k")),
    edu      = factor(edu,      levels = c("Less than highschool","Some highschool",
                                           "Highschool graduate","Some college",
                                           "College graduate or advanced")),
    homeown  = factor(homeown,  levels = c("Own","Rent","Others")),
    region   = factor(region,   levels = c("Northeast","Midwest","South","West"))
  )

########################################################################
## Save the per-person state file. Keep PWGTP raw over ALL adults so state
## population sums are complete (rho_a). Demographic columns are retained
## for cross-checks and may contain NA — do NOT na.omit this file.
########################################################################

dt_bystate <- dt %>%
  select(state, region, PWGTP, agegroup, sex, race_eth, income, edu, homeown)

## Sanity checks
message("Adult records: ", nrow(dt_bystate),
        " | states present: ", dplyr::n_distinct(dt_bystate$state))
message("Total adult population (sum PWGTP): ",
        format(round(sum(dt_bystate$PWGTP, na.rm = TRUE)), big.mark = ","))
if (any(is.na(dt_bystate$state)))
  warning(sum(is.na(dt_bystate$state)), " records have an unmapped state (ST) — check fips2usps.")
print(sort(table(dt_bystate$state)))

write.csv(dt_bystate, "../data/PUMS_2022_bystate.csv", row.names = FALSE)
message("Wrote ../data/PUMS_2022_bystate.csv")
