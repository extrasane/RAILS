#### VEHSS_data_eye.R — Download and assemble the VEHSS reference tables for
#### THREE eye phenotypes (AMD, glaucoma, cataract) from the CDC VEHSS portal
####
#### Multi-phenotype generalization of VEHSS_data.R (which stays as the AMD-
#### only script). For each phenotype it downloads, for every Age x Sex x
#### Race selection, the data the portal's "Export CSV" button hands out —
#### through the JSON endpoint the page itself calls,
####   https://ddt-vehss.cdc.gov/api/reportDataIds/<DataSource>/<Question>/?CountyFlag=N
####     &YearId=..&ResponseId=..&AgeId=..&GenderId=..&RaceId=..
####     &RiskFactorId=..&RiskFactorResponseId=..&DataValueTypeId=CRDPREV
#### — writes one csv per selection in the portal's export layout
#### (Export_<age>_<sex>_<race>.csv, one folder per phenotype) and merges them
#### into the tables VEHSS_comparison_eye.R reads:
####   VEHSS_<TAG>_overall.csv        40+ / both sexes / all races, by state (+ US)
####   VEHSS_<TAG>_combinations.csv   every other Age x Sex x Race selection
####   VEHSS_<TAG>_coverage.csv       the 4 x 3 x 5 grid, present / missing / derived
####
#### THE THREE INDICATORS (portal ids checked 2026-09-18 via /api/lookups and
#### /api/indicatorAllIds) and how their dimensions differ from AMD:
####
####   AMD       VEHSS Modeled Estimates (PREV) | QAMDM ~ R3_ALL "All AMD" |
####             YR11 = 2019 | RiskFactor RFPERS "All persons"
####             ages 40+ | 40-64 | 65-84 | 85+ ; sexes Both | Female | Male ;
####             races All | Black NH | Hispanic | White NH | Other (= all other
####             non-Hispanic groups, Asian included)
####
####   GLAUCOMA  VEHSS Modeled Estimates (PREV) | QGlauM ~ R5_ALL "All Glaucoma"
####             (open- or closed-angle glaucoma in either eye, exam-based, so
####             diagnosed + undiagnosed; glaucoma SUSPECTS not included) |
####             YR19 = 2022, the only modeled year | RiskFactor RFALL
####             ages ALL AGES | 0-39 | 40-64 | 65-84 | 85+  -> there is NO
####             "40 years and older" selection: the 40+ rows are DERIVED here
####             by summing cases and implied populations of the three 40+
####             bands (no confidence interval).
####             sexes and races as AMD.
####
####   CATARACT  Medicare Fee-for-Service AND Medicare Advantage claims
####             ("Medicare FFS MA") | QDXDC6 "Annual prevalence of diagnosed
####             cataracts" ~ R6_ALL "All Cataracts" (any Cat_6 subgroup:
####             senile, non-congenital, congenital, posterior capsular
####             opacity, pseudophakia, aphakia / other lens disorders) |
####             YR19 = 2022 (the portal id YR19 in the given URL is the year
####             2022, not 2019; the source has YR18 = 2021 and YR19 = 2022) |
####             RiskFactor RFALL "All patients" — the diabetes / hypertension
####             risk-factor splits (RFDM / RFHT) exist but are NOT used.
####             ages ALL | 0-17 | 18-39 | 18+ | 40-64 | 40+ | 65-84 | 65+ | 85+
####             (the four AMD codes all exist and are used);
####             sexes Both | Female | Male (same as AMD);
####             races All | White NH | Black NH | Hispanic | Asian | North
####             American Native | Other  -> the portal's "Other" EXCLUDES Asian
####             and AIAN, unlike the modeled sources and unlike AoU's "Other"
####             (NH Asian + Others). A comparable "Other" is therefore DERIVED
####             here as Asian + AIAN + Other (cases and implied populations
####             summed; no CI); the three portal groups are kept as well.
####             DENOMINATOR: Medicare beneficiaries, not the resident
####             population (implied 40+ population 52.8M vs ACS 162M). 65+ is
####             near-universal coverage; 40-64 is disability / ESRD
####             beneficiaries only, and MA encounter data are less complete
####             than FFS claims. Diagnosed prevalence, one code suffices.
####
#### Sample_Size in every source is the estimated NUMBER OF PEOPLE WITH THE
#### CONDITION in the selection (checked: state values sum to the US value and
#### the race groups sum to all races), so a selection's population is
#### Sample_Size / (Data_Value / 100).

suppressPackageStartupMessages({ library(tidyverse); library(httr); library(jsonlite) })

########################################################################
## Locations and options
########################################################################

LOCAL_VEHSS_DIR <- "C:/Work/2026 Spring/Sub RAILS/VEHSS"
VEHSS_DIR   <- if (dir.exists(dirname(LOCAL_VEHSS_DIR))) LOCAL_VEHSS_DIR else "VEHSS_exports"
EXTRA_FILES_AMD <- "C:/Work/2026 Spring/Export.csv"       # the original hand export (optional)
OUT_DIRS    <- unique(c(getwd(), "C:/Work/REpo/RAILS_Edit/Application/Subgroup RAILS", VEHSS_DIR))
OUT_DIRS    <- OUT_DIRS[vapply(OUT_DIRS, function(d) d == getwd() || d == VEHSS_DIR || dir.exists(d), logical(1))]

FORCE_DOWNLOAD <- FALSE     # TRUE = re-download every selection even if cached
DOWNLOAD       <- TRUE      # FALSE = only merge what is already on disk
API_PAUSE_SEC  <- 0.4       # politeness delay between API calls
API_BASE       <- "https://ddt-vehss.cdc.gov"

## Cataract reference: the source and year from the portal URL
## (DataSourceId=Medicare+FFS+MA, YearId=YR19). Alternatives with state x age
## x sex x race diagnosed-cataract data: "Medicare FFS" (2014-2022 incl. YR11 =
## 2019), "MEDICAID", "IRIS", "MSCANCC". The export folder name carries the
## source and year, so switching never mixes files from two sources.
CATARACT_SOURCE <- "Medicare FFS MA"
CATARACT_YEAR   <- "YR19"     # = 2022 in the portal's year lookup

## Race groups: portal code -> short code. The standard five are what
## VEHSS_comparison_eye.R uses; a phenotype may offer more (cataract), and
## `derive_races` lists derived groups (name = component groups) and
## `drop_races` the standalone portal levels removed once merged.
RACE_STD <- c(ALLRACE = "All", BLK = "Black", HISP = "Hispanic", WHT = "White", OTH = "Other")
RACE_MEDICARE <- c(ALLRACE = "All", BLK = "Black", HISP = "Hispanic", WHT = "White",
                   OTH = "Other (portal)", ASN = "Asian", AIAN = "AIAN")

PHENOTYPES <- list(
  amd = list(
    tag = "AMD", label = "Age-related macular degeneration (AMD)",
    data_source = "PREV", question = "QAMDM", response = "R3_ALL", year = "YR11",
    risk_factor = "RFPERS", rf_response = "RFTOT",
    response_pattern = "All Age-related macular degeneration",
    ages = c(AGE40PLUS = "40+", AGE4064 = "40-64", AGE6584 = "65-84", AGE85PLUS = "85+"),
    races = RACE_STD, derive_races = NULL, drop_races = NULL,
    export_dir = VEHSS_DIR, extra_files = EXTRA_FILES_AMD),
  glaucoma = list(
    tag = "GLAUCOMA", label = "Glaucoma",
    data_source = "PREV", question = "QGlauM", response = "R5_ALL", year = "YR19",
    risk_factor = "RFALL", rf_response = "RFTOT",
    response_pattern = "All Glaucoma",
    ages = c(AGE4064 = "40-64", AGE6584 = "65-84", AGE85PLUS = "85+"),   # no 40+ on the portal: derived
    races = RACE_STD, derive_races = NULL, drop_races = NULL,
    export_dir = file.path(VEHSS_DIR, "glaucoma"), extra_files = character(0)),
  cataract = list(
    tag = "CATARACT", label = paste0("Cataract (diagnosed, ", CATARACT_SOURCE, ")"),
    data_source = CATARACT_SOURCE, question = "QDXDC6", response = "R6_ALL", year = CATARACT_YEAR,
    risk_factor = "RFALL", rf_response = "RFTOT",
    response_pattern = "All Cataract",
    ages = c(AGE40PLUS = "40+", AGE4064 = "40-64", AGE6584 = "65-84", AGE85PLUS = "85+"),
    races = RACE_MEDICARE,
    ## derived race groups (name = components); the standard "Other" is the
    ## one the comparison scripts use (= AoU's Other: NH Asian + Others);
    ## "Other excl. Asian" folds North American Native into the portal's
    ## Other and is kept as an extra level together with "Asian"
    derive_races = list("Other" = c("Other (portal)", "Asian", "AIAN"),
                        "Other excl. Asian" = c("Other (portal)", "AIAN")),
    drop_races = c("AIAN", "Other (portal)"),
    export_dir = file.path(VEHSS_DIR, paste0("cataract_", gsub("[^A-Za-z0-9]+", "_", tolower(CATARACT_SOURCE)),
                                             "_", CATARACT_YEAR)),
    extra_files = character(0))
)

SEX_CODES  <- c(GALL = "Both", GF = "Female", GM = "Male")
AGE_LEVELS  <- c("40+", "40-64", "65-84", "85+")
SEX_LEVELS  <- unname(SEX_CODES); RACE_LEVELS <- unname(RACE_STD)

## Portal wording -> short codes (for reading exports by content)
AGE_MAP  <- c("40 years and older" = "40+", "40-64 years" = "40-64",
              "65-84 years" = "65-84", "85 years and older" = "85+")
SEX_MAP  <- c("Both sexes" = "Both", "Female" = "Female", "Male" = "Male")
race_map_of <- function(ph) {
  m <- c("All races" = "All", "Black, non-Hispanic" = "Black", "Hispanic, any race" = "Hispanic",
         "White, non-Hispanic" = "White", "Asian" = "Asian", "North American Native" = "AIAN",
         "Other" = if ("Other (portal)" %in% ph$races) "Other (portal)" else "Other")
  m
}

## File-name slugs for the race short codes (and back)
RACE_SLUG <- c("All" = "all", "Black" = "black", "Hispanic" = "hispanic", "White" = "white",
               "Other" = "other", "Other (portal)" = "otherportal", "Asian" = "asian", "AIAN" = "aian")
file_slug <- function(age, sex, race)
  paste0("Export_", sub("\\+", "plus", age), "_", tolower(sex), "_", RACE_SLUG[race], ".csv")

########################################################################
## Portal API helpers
########################################################################

ua <- user_agent("R / RAILS VEHSS comparison (research use)")
## The claims-source endpoints (Medicare) answer slowly and intermittently
## with HTTP 503: retry up to 8 times with exponential back-off (2, 4, 8, ...
## seconds, capped at 60) and allow 120 s per request.
api_get <- function(path, query = NULL) {
  r <- RETRY("GET", paste0(API_BASE, path), query = query, ua, timeout(120),
             times = 8, pause_base = 2, pause_cap = 60, pause_min = 1, quiet = TRUE)
  if (http_error(r)) stop("VEHSS API error ", status_code(r), " for ", path)
  txt <- content(r, as = "text", encoding = "UTF-8")
  if (!startsWith(trimws(txt), "[") && !startsWith(trimws(txt), "{"))
    stop("VEHSS API did not return JSON for ", path, ": ", substr(txt, 1, 200))
  fromJSON(txt)
}
enc <- function(x) URLencode(x, reserved = TRUE)     # data-source ids may contain spaces

message("Fetching portal lookups ...")
lookups   <- api_get("/api/lookups")
footnotes <- api_get("/api/footnotes")
lk <- function(type, id) {
  m <- lookups[lookups$type == type, ]
  out <- m$name[match(id, m$id)]
  ifelse(is.na(out), id, out)
}
loc_tab <- lookups[lookups$type == "Location", c("id", "name", "abbr")]
fn_text <- function(fs) {
  if (is.na(fs) || !nzchar(fs)) return("")
  hit <- footnotes$ft[footnotes$fs == fs]
  if (length(hit)) hit[1] else fs
}

## Which dimension ids the portal actually has for a phenotype
available_dims <- function(ph) {
  d <- tryCatch(api_get(paste0("/api/indicatorAllIds/", enc(ph$data_source), "/",
                               ph$question, "~", ph$response)),
                error = function(e) NULL)
  if (is.null(d) || length(d) == 0) return(NULL)
  as.data.frame(d)
}

## One selection -> data frame in the portal's export layout
fetch_selection <- function(ph, age_id, sex_id, race_id) {
  q <- list(CountyFlag = "N", YearId = ph$year, ResponseId = ph$response,
            AgeId = age_id, GenderId = sex_id, RaceId = race_id,
            RiskFactorId = ph$risk_factor, RiskFactorResponseId = ph$rf_response,
            DataValueTypeId = "CRDPREV")
  d <- api_get(paste0("/api/reportDataIds/", enc(ph$data_source), "/", ph$question, "/"), query = q)
  if (length(d) == 0 || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  li <- match(d$loc, loc_tab$id)
  tibble(
    YearStart = lk("Year", d$yr), YearEnd = lk("Year", d$yr),
    LocationAbbr = loc_tab$abbr[li], LocationDesc = loc_tab$name[li], CountyName = "",
    DataSource = ph$data_source, Topic = "", Category = "",
    Question   = lk("Question", ph$question),
    Response   = lk("Response", d$rs),
    Age        = lk("AgeGroup", d$ag),
    Sex        = lk("Sex", d$ge),
    Race_Ethnicity = lk("RaceEthnicity", d$re),
    Risk_Factor    = lk("RiskFactor", d$rf),
    Risk_Factor_Repsonse = lk("RiskFactorResponse", d$rfr),   # (sic) the portal's spelling
    Data_Value_Unit = d$dvu, Data_Value_Type = lk("DataValueType", "CRDPREV"),
    Data_Value = d$dv, Data_Value_Footnote_Symbol = d$fs,
    Data_Value_Footnote = vapply(d$fs, fn_text, character(1)),
    Low_Confidence_Limit = d$lci, High_Confidence_Limit = d$hci, Sample_Size = d$ss,
    FootnoteText = "", Url = "", FootnoteType = ""
  ) %>% arrange(LocationAbbr != "US", LocationAbbr)
}

########################################################################
## Reading exports by CONTENT (never trusting the file name)
########################################################################

harmonize <- function(x, map, what) {
  out <- unname(map[x])
  bad <- unique(x[is.na(out) & !is.na(x)])
  if (length(bad)) {
    warning("Unrecognized ", what, " label(s) kept verbatim: ", paste(bad, collapse = " | "))
    out[is.na(out)] <- x[is.na(out)]
  }
  out
}
claim_from_name <- function(f) {
  b <- tolower(tools::file_path_sans_ext(basename(f)))
  if (b == "export") return(c(age = "40+", sex = "Both", race = "All"))
  m <- regmatches(b, regexec("^export_([0-9]+-[0-9]+|[0-9]+\\+|[0-9]+plus)_(both|female|male|all)_([a-z]+)$", b))[[1]]
  if (length(m) == 0) return(c(age = NA, sex = NA, race = NA))
  inv <- setNames(names(RACE_SLUG), RACE_SLUG); inv["others"] <- "Other"
  c(age = sub("plus$", "+", m[2]),
    sex = unname(c(both = "Both", all = "Both", female = "Female", male = "Male")[m[3]]),
    race = unname(inv[m[4]]))
}
read_export <- function(f, ph) {
  d <- readr::read_csv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c")) %>%
    filter(!is.na(LocationAbbr), nzchar(LocationAbbr), !is.na(Age))
  if (nrow(d) == 0) { warning("No data rows in ", basename(f)); return(NULL) }
  ok_ind <- grepl(ph$response_pattern, d$Response, ignore.case = TRUE) & grepl("Crude", d$Data_Value_Type)
  if (!all(ok_ind))
    warning(basename(f), ": ", sum(!ok_ind), " rows are not '", ph$response_pattern, " / crude prevalence' — dropped")
  d <- d[ok_ind, ]
  d %>%
    transmute(source_file = basename(f), file_mtime = file.info(f)$mtime,
              data_source = DataSource, year = YearStart,
              age  = harmonize(Age, AGE_MAP, "age"),
              sex  = harmonize(Sex, SEX_MAP, "sex"),
              race = harmonize(Race_Ethnicity, race_map_of(ph), "race"),
              state = LocationAbbr, state_name = LocationDesc,
              prev_pct = suppressWarnings(as.numeric(Data_Value)),
              lb_pct   = suppressWarnings(as.numeric(Low_Confidence_Limit)),
              ub_pct   = suppressWarnings(as.numeric(High_Confidence_Limit)),
              cases    = suppressWarnings(as.numeric(Sample_Size)),
              footnote = Data_Value_Footnote)
}
## A cached file is reusable only if it holds exactly the selection its name
## claims AND comes from this phenotype's source and year (a stale export from
## another source / year is re-downloaded)
file_matches_name <- function(f, ph) {
  d <- tryCatch(suppressWarnings(read_export(f, ph)), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(FALSE)
  cl  <- claim_from_name(f)
  got <- unique(paste(d$age, d$sex, d$race))
  length(got) == 1 && identical(got, paste(cl["age"], cl["sex"], cl["race"])) &&
    all(d$year == lk("Year", ph$year)) &&
    all(d$data_source %in% c(ph$data_source, lk("DataSource", ph$data_source)))
}

########################################################################
## Per-phenotype: download, merge, derive, write tables
########################################################################

dim_notes <- list()

process_phenotype <- function(ph) {
  message("\n==================== ", ph$label, " ====================")
  dims <- available_dims(ph)
  if (is.null(dims)) {
    warning(ph$label, ": the portal has no data for ", ph$data_source, " / ", ph$question, "~", ph$response,
            " — skipped.", immediate. = TRUE)
    return(invisible(NULL))
  }
  dims <- dims[dims$dv == "CRDPREV", ]
  if (!ph$year %in% dims$yr)
    stop(ph$label, ": year ", ph$year, " not offered; available: ", paste(sort(unique(dims$yr)), collapse = ", "))
  dims_y <- dims[dims$yr == ph$year & dims$rs == ph$response, ]
  if (!ph$risk_factor %in% dims_y$rf) {
    old <- ph$risk_factor; ph$risk_factor <- sort(unique(dims_y$rf))[1]
    message("  RiskFactorId ", old, " not offered; using ", ph$risk_factor)
  }
  dims_y <- dims_y[dims_y$rf == ph$risk_factor, ]
  offered <- list(ages = lk("AgeGroup", sort(unique(dims_y$ag))), sexes = lk("Sex", sort(unique(dims_y$ge))),
                  races = lk("RaceEthnicity", sort(unique(dims_y$re))), years = lk("Year", sort(unique(dims$yr))),
                  risk_factors = lk("RiskFactor", sort(unique(dims$rf))))
  ages_ok <- names(ph$ages)[names(ph$ages) %in% dims_y$ag]
  if (length(ages_ok) < length(ph$ages))
    message("  age codes not offered for ", lk("Year", ph$year), ": ",
            paste(setdiff(names(ph$ages), ages_ok), collapse = ", "), " (derived below if 40+)")
  ph$ages <- ph$ages[ages_ok]
  races_ok <- names(ph$races)[names(ph$races) %in% dims_y$re]
  ph$races <- ph$races[races_ok]
  message("  source ", ph$data_source, " | ", ph$question, "~", ph$response, " | year ", lk("Year", ph$year),
          " | RF ", ph$risk_factor, "\n  portal offers: years ", paste(offered$years, collapse = ", "),
          "\n                 ages ", paste(offered$ages, collapse = " | "),
          "\n                 sexes ", paste(offered$sexes, collapse = " | "),
          "\n                 races ", paste(offered$races, collapse = " | "),
          "\n                 risk factors ", paste(offered$risk_factors, collapse = " | "),
          "\n  used here:     ages ", paste(ph$ages, collapse = ", "), if (!"40+" %in% ph$ages) " (+ 40+ derived)" else "",
          " | races ", paste(ph$races, collapse = ", "),
          if (!is.null(ph$derive_races)) paste0(" (+ derived: ", paste(vapply(names(ph$derive_races), function(n)
            paste0(n, " = ", paste(ph$derive_races[[n]], collapse = " + ")), character(1)), collapse = "; "),
            "; dropped: ", paste(ph$drop_races, collapse = ", "), ")") else "")
  dim_notes[[ph$tag]] <<- c(list(label = ph$label, source = ph$data_source, year = lk("Year", ph$year),
                                 used_ages = paste(ph$ages, collapse = ", "), used_races = paste(ph$races, collapse = ", ")),
                            offered)

  dir.create(ph$export_dir, showWarnings = FALSE, recursive = TRUE)
  grid <- expand_grid(age_id = names(ph$ages), sex_id = names(SEX_CODES), race_id = names(ph$races)) %>%
    mutate(age = ph$ages[age_id], sex = SEX_CODES[sex_id], race = ph$races[race_id],
           file = file.path(ph$export_dir, file_slug(age, sex, race)))

  ## ---- 1. download ------------------------------------------------------
  if (DOWNLOAD) {
    message("  downloading ", nrow(grid), " Age x Sex x Race selections into ", ph$export_dir, " ...")
    status <- character(nrow(grid))
    for (i in seq_len(nrow(grid))) {
      g <- grid[i, ]
      tag <- paste(g$age, g$sex, g$race, sep = " / ")
      if (!FORCE_DOWNLOAD && file.exists(g$file)) {
        if (file_matches_name(g$file, ph)) { status[i] <- "cached"; next }
        warning("Cached ", basename(g$file), " does not contain ", tag, " from ", ph$data_source, " ",
                lk("Year", ph$year), " — re-downloading.", immediate. = TRUE)
      }
      d <- tryCatch(fetch_selection(ph, g$age_id, g$sex_id, g$race_id),
                    error = function(e) { warning(tag, ": ", conditionMessage(e), immediate. = TRUE); "error" })
      if (identical(d, "error")) { status[i] <- "ERROR"; next }
      if (is.null(d)) { status[i] <- "no data"; next }
      write_excel_csv(d, g$file, na = "")
      status[i] <- "downloaded"
      Sys.sleep(API_PAUSE_SEC)
    }
    message("  download summary: ", paste(sprintf("%s = %d", names(table(status)), table(status)), collapse = ", "))
  }

  ## ---- 2. read everything in export_dir (+ extra files) by content -------
  files <- c(list.files(ph$export_dir, pattern = "^Export.*[.]csv$", full.names = TRUE), ph$extra_files)
  files <- files[file.exists(files)]
  if (length(files) == 0) { warning(ph$label, ": no exports found in ", ph$export_dir); return(invisible(NULL)) }
  raw <- bind_rows(lapply(files, read_export, ph = ph))
  ## keep only this phenotype's source + year (stale files from another source are ignored)
  keep_src <- raw$year == lk("Year", ph$year) & raw$data_source %in% c(ph$data_source, lk("DataSource", ph$data_source))
  if (any(!keep_src))
    message("  ", sum(!keep_src), " rows from other sources / years ignored (",
            paste(unique(raw$source_file[!keep_src]), collapse = ", "), ")")
  raw <- raw[keep_src, ]

  chk <- raw %>% group_by(source_file) %>%
    summarise(content = paste(unique(paste(age, sex, race, sep = " / ")), collapse = " ; "),
              n_sel = n_distinct(paste(age, sex, race)), .groups = "drop") %>%
    rowwise() %>%
    mutate(claimed = { cl <- claim_from_name(source_file)
                       if (all(is.na(cl))) NA_character_ else paste(cl["age"], cl["sex"], cl["race"], sep = " / ") },
           status = case_when(n_sel > 1 ~ "MULTIPLE", is.na(claimed) ~ "unparsed name",
                              claimed == content ~ "ok", TRUE ~ "NAME DISAGREES WITH CONTENT")) %>% ungroup()
  bad <- chk %>% filter(status != "ok")
  if (nrow(bad)) { message("  files whose name and content differ (content is used):"); print(bad, n = Inf, width = Inf) }

  dedup <- raw %>% group_by(age, sex, race) %>%
    group_modify(function(g, key) {
      fs <- unique(g$source_file); if (length(fs) == 1) return(g)
      named <- fs[vapply(fs, function(f) { cl <- claim_from_name(f)
        !any(is.na(cl)) && identical(unname(cl), c(key$age, key$sex, key$race)) }, logical(1))]
      keep <- if (length(named)) named[1] else
        g %>% group_by(source_file) %>% summarise(m = max(file_mtime), .groups = "drop") %>%
          slice_max(m, n = 1, with_ties = FALSE) %>% pull(source_file)
      g %>% filter(source_file == keep)
    }) %>% ungroup()

  tab <- dedup %>%
    mutate(prev = prev_pct / 100, lb = lb_pct / 100, ub = ub_pct / 100, pop = cases / prev,
           derived = FALSE)

  ## ---- 3a. derive race groups when the portal splits "Other" ---------------
  ## Each derived group is the sum of its component groups (cases and implied
  ## populations). When a component is suppressed (AIAN often is, in small
  ## states) the same quantity is recovered as All races minus the groups
  ## OUTSIDE the derived one — exact up to the portal's rounding, since the
  ## race groups add up to all races. No confidence interval either way.
  ## The standalone AIAN and portal-"Other" levels are then dropped.
  if (!is.null(ph$derive_races)) {
    all_groups <- setdiff(unique(as.character(tab$race)), "All")
    for (new_name in names(ph$derive_races)) {
      comps <- ph$derive_races[[new_name]]
      minus <- setdiff(all_groups, comps)
      d_sum <- tab %>% filter(race %in% comps) %>%
        group_by(age, sex, state, state_name, year) %>%
        summarise(ok_sum = n_distinct(race) == length(comps) & !any(is.na(prev)),
                  cases_sum = sum(cases), pop_sum = sum(pop), .groups = "drop")
      d_res <- tab %>% filter(race %in% c("All", minus)) %>%
        group_by(age, sex, state, state_name, year) %>%
        summarise(ok_res = n_distinct(race) == length(minus) + 1 & !any(is.na(prev)),
                  cases_res = sum(cases[race == "All"]) - sum(cases[race %in% minus]),
                  pop_res   = sum(pop[race == "All"])   - sum(pop[race %in% minus]), .groups = "drop")
      d_new <- full_join(d_sum, d_res, by = c("age", "sex", "state", "state_name", "year")) %>%
        mutate(ok_sum = coalesce(ok_sum, FALSE), ok_res = coalesce(ok_res, FALSE),
               race = new_name,
               cases = case_when(ok_sum ~ cases_sum, ok_res ~ pmax(cases_res, 0), TRUE ~ NA_real_),
               pop   = case_when(ok_sum ~ pop_sum,   ok_res ~ pmax(pop_res, 0),   TRUE ~ NA_real_),
               prev = cases / pop, lb = NA_real_, ub = NA_real_,
               prev_pct = 100 * prev, lb_pct = NA_real_, ub_pct = NA_real_,
               footnote = case_when(
                 ok_sum ~ paste0("Derived: ", new_name, " = ", paste(comps, collapse = " + "),
                                 " (cases and implied populations summed; no CI)"),
                 ok_res ~ paste0("Derived: ", new_name, " = All races - ", paste(minus, collapse = " - "),
                                 " (a component group is suppressed; no CI)"),
                 TRUE   ~ paste0("Derived ", new_name, ": components and the residual both unavailable")),
               source_file = "derived from race groups", derived = TRUE, data_source = ph$data_source) %>%
        select(age, sex, state, state_name, year, race, cases, pop, prev, lb, ub, prev_pct, lb_pct, ub_pct,
               footnote, source_file, derived, data_source, ok_sum, ok_res)
      message("  '", new_name, "' race rows DERIVED for ", n_distinct(paste(d_new$age, d_new$sex)),
              " age x sex selections: ", sum(d_new$ok_sum), " state values as ", paste(comps, collapse = " + "),
              ", ", sum(!d_new$ok_sum & d_new$ok_res), " as All - ", paste(minus, collapse = " - "),
              ", ", sum(is.na(d_new$prev)), " missing")
      tab <- bind_rows(tab, d_new %>% select(-ok_sum, -ok_res))
    }
    if (length(ph$drop_races)) {
      message("  standalone race levels dropped (merged into the derived groups): ", paste(ph$drop_races, collapse = ", "))
      tab <- tab %>% filter(!race %in% ph$drop_races)
    }
  }

  ## ---- 3b. derive the 40+ rows from the bands when the portal has none ----
  if (!"40+" %in% tab$age) {
    bands <- c("40-64", "65-84", "85+")
    d40 <- tab %>% filter(age %in% bands) %>%
      group_by(sex, race, state, state_name) %>%
      summarise(n_bands = n_distinct(age), any_na = any(is.na(prev)),
                cases = sum(cases), pop = sum(pop), year = paste(sort(unique(year)), collapse = "/"),
                source_file = "derived from age bands", .groups = "drop") %>%
      mutate(age = "40+",
             cases = ifelse(n_bands == length(bands) & !any_na, cases, NA_real_),
             pop   = ifelse(n_bands == length(bands) & !any_na, pop,   NA_real_),
             prev  = cases / pop, lb = NA_real_, ub = NA_real_,
             prev_pct = 100 * prev, lb_pct = NA_real_, ub_pct = NA_real_,
             footnote = ifelse(is.na(prev), "Derived 40+: a band is missing or suppressed",
                               "Derived: 40+ aggregated from 40-64 / 65-84 / 85+ (cases and implied population summed; no CI)"),
             derived = TRUE, data_source = ph$data_source) %>%
      select(-n_bands, -any_na)
    message("  40+ rows DERIVED from the age bands for ", n_distinct(paste(d40$sex, d40$race)),
            " sex x race selections (", sum(is.na(d40$prev)), " state values missing)")
    tab <- bind_rows(tab, d40)
  }

  tab <- tab %>%
    mutate(age  = factor(age,  levels = union(AGE_LEVELS,  unique(age))),
           sex  = factor(sex,  levels = union(SEX_LEVELS,  unique(sex))),
           race = factor(race, levels = union(RACE_LEVELS, unique(race)))) %>%
    select(age, sex, race, state, state_name, year, data_source, prev, lb, ub, prev_pct, lb_pct, ub_pct,
           cases, pop, derived, footnote, source_file) %>%
    arrange(age, sex, race, state != "US", state)

  is_overall <- tab$age == "40+" & tab$sex == "Both" & tab$race == "All"
  overall <- tab[is_overall, ]; combos <- tab[!is_overall, ]
  if (nrow(overall) == 0) warning(ph$label, ": the overall (40+ / Both / All) selection is missing.")

  coverage <- expand_grid(age = AGE_LEVELS, sex = SEX_LEVELS, race = RACE_LEVELS) %>%
    left_join(tab %>% group_by(age, sex, race) %>%
                summarise(n_states = sum(state != "US"), n_suppressed = sum(is.na(prev) & state != "US"),
                          has_US = any(state == "US"), derived = any(derived), .groups = "drop") %>%
                mutate(across(c(age, sex, race), as.character)), by = c("age", "sex", "race")) %>%
    mutate(present = !is.na(n_states)) %>% arrange(desc(present), age, sex, race)
  message("  coverage: ", sum(coverage$present), " of ", nrow(coverage), " standard selections present (",
          sum(coverage$derived, na.rm = TRUE), " derived); suppressed state values: ",
          sum(coverage$n_suppressed, na.rm = TRUE),
          "; extra race groups kept: ", paste(setdiff(unique(as.character(tab$race)), RACE_LEVELS), collapse = ", "))
  us <- overall %>% filter(state == "US")
  if (nrow(us)) message("  national ", ph$label, " 40+: ", sprintf("%.2f%%", 100 * us$prev),
                        if (!is.na(us$lb)) sprintf(" (%.2f-%.2f)", 100 * us$lb, 100 * us$ub) else " (derived, no CI)",
                        " | implied 40+ population ", format(round(us$pop), big.mark = ","))

  for (d in OUT_DIRS) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    write_excel_csv(overall,  file.path(d, paste0("VEHSS_", ph$tag, "_overall.csv")))
    write_excel_csv(combos,   file.path(d, paste0("VEHSS_", ph$tag, "_combinations.csv")))
    write_excel_csv(coverage, file.path(d, paste0("VEHSS_", ph$tag, "_coverage.csv")))
  }
  message("  wrote VEHSS_", ph$tag, "_{overall,combinations,coverage}.csv -> ", paste(OUT_DIRS, collapse = " ; "))
  invisible(tab)
}

results <- lapply(PHENOTYPES, process_phenotype)

########################################################################
## Dimension differences between the three indicators (also printed above)
########################################################################

message("\n==================== DIMENSION DIFFERENCES vs AMD ====================")
for (tg in names(dim_notes)) {
  n <- dim_notes[[tg]]
  message(sprintf("%-9s %s | %s | ages offered: %s | ages used: %s | sexes: %s | races offered: %s",
                  tg, n$source, n$year, paste(n$ages, collapse = " / "), n$used_ages,
                  paste(n$sexes, collapse = " / "), paste(n$races, collapse = " / ")))
}
message("Sexes are identical across the three (Both / Female / Male). Age bands 40-64 / 65-84 / 85+ exist in all ",
        "three; '40 years and older' is missing for glaucoma (derived) and present for AMD and cataract. Races: the ",
        "Medicare cataract source splits Asian and North American Native out of 'Other' (a comparable 'Other' is ",
        "derived); the cataract denominator is Medicare beneficiaries, not the resident population.")
message("\nDone: ", paste(names(results)[!vapply(results, is.null, logical(1))], collapse = ", "))
