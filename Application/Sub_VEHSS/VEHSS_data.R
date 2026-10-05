#### VEHSS_data.R — Download and assemble the CDC VEHSS AMD reference tables
####
#### Downloads, for every Age x Sex x Race/ethnicity selection, the same data
#### the VEHSS Data Portal (https://ddt-vehss.cdc.gov) hands out through its
#### "Export CSV" button, writes one csv per selection in the portal's export
#### layout (Export_<age>_<sex>_<race>.csv in VEHSS_DIR, exactly like the files
#### exported by hand), and merges them into the two tables VEHSS_comparison.R
#### reads.
####
#### Why download instead of exporting by hand: the portal's Export button
#### builds the csv from whatever data the page last fetched, so changing a
#### dropdown and exporting can hand back the PREVIOUS selection (three of the
#### first eight hand-made exports carried a different Age than their name).
#### The page itself fetches its numbers from a JSON endpoint,
####   https://ddt-vehss.cdc.gov/api/reportDataIds/PREV/QAMDM/?CountyFlag=N
####     &YearId=YR11&ResponseId=R3_ALL&AgeId=<age>&GenderId=<sex>&RaceId=<race>
####     &RiskFactorId=RFPERS&RiskFactorResponseId=RFTOT&DataValueTypeId=CRDPREV
#### which this script calls directly with the codes the portal's dropdowns use:
####   AgeId    AGE40PLUS | AGE4064 | AGE6584 | AGE85PLUS
####   GenderId GALL (both sexes) | GF | GM
####   RaceId   ALLRACE | BLK (Black, non-Hispanic) | HISP (Hispanic, any race)
####            | WHT (White, non-Hispanic) | OTH (Other)
#### Labels (age / sex / race / location / year / footnotes) come from the
#### portal's /api/lookups and /api/footnotes, so the csvs match the manual
#### exports column for column.
####
#### Indicator: "Prevalence of AMD" (QAMDM), response "All Age-related macular
#### degeneration (AMD)" (R3_ALL), VEHSS Modeled Estimates (PREV), crude
#### prevalence (%), year YR11 = 2019, state + national (CountyFlag = N).
#### Sample_Size is the estimated NUMBER OF PEOPLE WITH AMD in the selection
#### (state values sum to the US value), so its population is
#### Sample_Size / (Data_Value / 100). Suppressed values arrive as NA with a
#### footnote symbol.
####
#### Steps
####   1. download every selection not yet cached (FORCE_DOWNLOAD re-downloads
####      all; a cached file whose CONTENT does not match its name — a stale
####      hand export — is re-downloaded and overwritten, with a warning),
####   2. read all Export*.csv in VEHSS_DIR (+ EXTRA_FILES), classify by content,
####      de-duplicate, harmonize labels,
####   3. write to OUT_DIRS:
####        VEHSS_AMD_overall.csv        40+ / both sexes / all races, by state (+ US)
####        VEHSS_AMD_combinations.csv   every other Age x Sex x Race selection
####        VEHSS_AMD_coverage.csv       the 4 x 3 x 5 grid, present / missing

suppressPackageStartupMessages({ library(tidyverse); library(httr); library(jsonlite) })

########################################################################
## Locations and options
########################################################################

## On the local machine the exports live in the 2026 Spring folder; anywhere
## else (e.g. sourced from VEHSS_comparison.R inside the AoU Workbench) they
## go to ./VEHSS_exports. The merged tables are always written to the working
## directory (where VEHSS_comparison.R looks for them) and to VEHSS_DIR.
LOCAL_VEHSS_DIR <- "C:/Work/2026 Spring/Sub RAILS/VEHSS"
VEHSS_DIR   <- if (dir.exists(dirname(LOCAL_VEHSS_DIR))) LOCAL_VEHSS_DIR else "VEHSS_exports"
EXTRA_FILES <- "C:/Work/2026 Spring/Export.csv"             # the original hand export (optional)
OUT_DIRS    <- unique(c(getwd(), "C:/Work/REpo/RAILS_Edit/Application/Subgroup RAILS", VEHSS_DIR))
OUT_DIRS    <- OUT_DIRS[vapply(OUT_DIRS, function(d) d == getwd() || d == VEHSS_DIR || dir.exists(d), logical(1))]

FORCE_DOWNLOAD <- FALSE     # TRUE = re-download every selection even if cached
DOWNLOAD       <- TRUE      # FALSE = only merge what is already in VEHSS_DIR
API_PAUSE_SEC  <- 0.4       # politeness delay between API calls

API_BASE  <- "https://ddt-vehss.cdc.gov"
API_QUERY <- list(CountyFlag = "N", YearId = "YR11", ResponseId = "R3_ALL",
                  RiskFactorId = "RFPERS", RiskFactorResponseId = "RFTOT",
                  DataValueTypeId = "CRDPREV")
DATA_SOURCE_ID <- "PREV"; QUESTION_ID <- "QAMDM"

## Portal codes -> short codes used in file names and in VEHSS_comparison.R
AGE_CODES  <- c(AGE40PLUS = "40+", AGE4064 = "40-64", AGE6584 = "65-84", AGE85PLUS = "85+")
SEX_CODES  <- c(GALL = "Both", GF = "Female", GM = "Male")
RACE_CODES <- c(ALLRACE = "All", BLK = "Black", HISP = "Hispanic", WHT = "White", OTH = "Other")
AGE_LEVELS  <- unname(AGE_CODES); SEX_LEVELS <- unname(SEX_CODES); RACE_LEVELS <- unname(RACE_CODES)

## Portal wording -> short codes (for reading hand-made exports by content)
AGE_MAP  <- c("40 years and older" = "40+", "40-64 years" = "40-64",
              "65-84 years" = "65-84", "85 years and older" = "85+")
SEX_MAP  <- c("Both sexes" = "Both", "Female" = "Female", "Male" = "Male")
RACE_MAP <- c("All races" = "All", "Black, non-Hispanic" = "Black",
              "Hispanic, any race" = "Hispanic", "White, non-Hispanic" = "White", "Other" = "Other")

## File name for a selection: Export_40plus_both_all.csv, Export_40-64_female_black.csv, ...
file_slug <- function(age, sex, race)
  paste0("Export_", sub("\\+", "plus", age), "_", tolower(sex), "_", tolower(race), ".csv")

########################################################################
## Portal API helpers
########################################################################

ua <- user_agent("R / RAILS VEHSS comparison (research use)")

api_get <- function(path, query = NULL) {
  r <- RETRY("GET", paste0(API_BASE, path), query = query, ua, times = 3, pause_base = 2, quiet = TRUE)
  if (http_error(r)) stop("VEHSS API error ", status_code(r), " for ", path)
  txt <- content(r, as = "text", encoding = "UTF-8")
  if (!startsWith(trimws(txt), "[") && !startsWith(trimws(txt), "{"))
    stop("VEHSS API did not return JSON for ", path, ": ", substr(txt, 1, 200))
  fromJSON(txt)
}

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
  if (!nzchar(fs)) return("")
  hit <- footnotes$ft[footnotes$fs == fs]
  if (length(hit)) hit[1] else fs
}
year_of <- function(yr) lk("Year", yr)

## One selection -> data frame in the portal's export layout
fetch_selection <- function(age_id, sex_id, race_id) {
  q <- c(API_QUERY, list(AgeId = age_id, GenderId = sex_id, RaceId = race_id))
  d <- api_get(paste0("/api/reportDataIds/", DATA_SOURCE_ID, "/", QUESTION_ID, "/"), query = q)
  if (length(d) == 0 || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  li <- match(d$loc, loc_tab$id)
  tibble(
    YearStart  = year_of(d$yr), YearEnd = year_of(d$yr),
    LocationAbbr = loc_tab$abbr[li], LocationDesc = loc_tab$name[li], CountyName = "",
    DataSource = DATA_SOURCE_ID, Topic = "", Category = "",
    Question   = lk("Question", QUESTION_ID),
    Response   = lk("Response", d$rs),
    Age        = lk("AgeGroup", d$ag),
    Sex        = lk("Sex", d$ge),
    Race_Ethnicity = lk("RaceEthnicity", d$re),
    Risk_Factor    = lk("RiskFactor", d$rf),
    Risk_Factor_Repsonse = lk("RiskFactorResponse", d$rfr),   # (sic) the portal's spelling
    Data_Value_Unit = d$dvu, Data_Value_Type = lk("DataValueType", API_QUERY$DataValueTypeId),
    Data_Value = d$dv, Data_Value_Footnote_Symbol = d$fs,
    Data_Value_Footnote = vapply(d$fs, fn_text, character(1)),
    Low_Confidence_Limit = d$lci, High_Confidence_Limit = d$hci, Sample_Size = d$ss,
    FootnoteText = "", Url = "", FootnoteType = ""
  ) %>% arrange(LocationAbbr != "US", LocationAbbr)
}

########################################################################
## Reading exports (downloaded or hand-made) by CONTENT
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

## Selection claimed by a file name (NA when the name has no such pattern)
claim_from_name <- function(f) {
  b <- tolower(tools::file_path_sans_ext(basename(f)))
  if (b == "export") return(c(age = "40+", sex = "Both", race = "All"))
  m <- regmatches(b, regexec("^export_([0-9]+-[0-9]+|[0-9]+\\+|[0-9]+plus)_(both|female|male|all)_([a-z]+)$", b))[[1]]
  if (length(m) == 0) return(c(age = NA, sex = NA, race = NA))
  age  <- sub("plus$", "+", m[2])
  sex  <- c(both = "Both", all = "Both", female = "Female", male = "Male")[m[3]]
  race <- c(all = "All", black = "Black", hispanic = "Hispanic", white = "White",
            other = "Other", others = "Other")[m[4]]
  c(age = age, sex = unname(sex), race = unname(race))
}

read_export <- function(f) {
  d <- readr::read_csv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c")) %>%
    filter(!is.na(LocationAbbr), nzchar(LocationAbbr), !is.na(Age))   # drops the footnote rows
  if (nrow(d) == 0) { warning("No data rows in ", basename(f)); return(NULL) }
  ok_ind <- grepl("All Age-related macular degeneration", d$Response) & grepl("Crude", d$Data_Value_Type)
  if (!all(ok_ind))
    warning(basename(f), ": ", sum(!ok_ind), " rows are not 'All AMD / crude prevalence' — dropped")
  d <- d[ok_ind, ]
  d %>%
    transmute(source_file = basename(f), file_mtime = file.info(f)$mtime,
              year = as.integer(YearStart),
              age  = harmonize(Age, AGE_MAP, "age"),
              sex  = harmonize(Sex, SEX_MAP, "sex"),
              race = harmonize(Race_Ethnicity, RACE_MAP, "race"),
              state = LocationAbbr, state_name = LocationDesc,
              prev_pct = suppressWarnings(as.numeric(Data_Value)),
              lb_pct   = suppressWarnings(as.numeric(Low_Confidence_Limit)),
              ub_pct   = suppressWarnings(as.numeric(High_Confidence_Limit)),
              cases    = suppressWarnings(as.numeric(Sample_Size)),
              footnote = Data_Value_Footnote)
}

## Does a cached file contain exactly the selection its name claims?
file_matches_name <- function(f) {
  d <- tryCatch(suppressWarnings(read_export(f)), error = function(e) NULL)
  if (is.null(d)) return(FALSE)
  cl  <- claim_from_name(f)
  got <- unique(paste(d$age, d$sex, d$race))
  length(got) == 1 && identical(got, paste(cl["age"], cl["sex"], cl["race"]))
}

########################################################################
## 1. Download every selection into VEHSS_DIR
########################################################################

dir.create(VEHSS_DIR, showWarnings = FALSE, recursive = TRUE)
grid <- expand_grid(age_id = names(AGE_CODES), sex_id = names(SEX_CODES), race_id = names(RACE_CODES)) %>%
  mutate(age = AGE_CODES[age_id], sex = SEX_CODES[sex_id], race = RACE_CODES[race_id],
         file = file.path(VEHSS_DIR, file_slug(age, sex, race)))

if (DOWNLOAD) {
  message("Downloading ", nrow(grid), " Age x Sex x Race selections from the VEHSS portal ...")
  dl_log <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    tag <- paste(g$age, g$sex, g$race, sep = " / ")
    if (!FORCE_DOWNLOAD && file.exists(g$file)) {
      if (file_matches_name(g$file)) { dl_log[[i]] <- tibble(selection = tag, status = "cached"); next }
      warning("Cached ", basename(g$file), " does not contain ", tag,
              " (stale hand export) — re-downloading and overwriting.", immediate. = TRUE)
    }
    d <- tryCatch(fetch_selection(g$age_id, g$sex_id, g$race_id),
                  error = function(e) { warning(tag, ": ", conditionMessage(e), immediate. = TRUE); "error" })
    if (identical(d, "error")) { dl_log[[i]] <- tibble(selection = tag, status = "ERROR"); next }
    if (is.null(d)) {
      message("  [", i, "/", nrow(grid), "] ", tag, ": no data on the portal")
      dl_log[[i]] <- tibble(selection = tag, status = "no data"); next
    }
    write_excel_csv(d, g$file, na = "")
    message("  [", i, "/", nrow(grid), "] ", tag, ": ", nrow(d), " rows -> ", basename(g$file))
    dl_log[[i]] <- tibble(selection = tag, status = paste0("downloaded (", nrow(d), " rows)"))
    Sys.sleep(API_PAUSE_SEC)
  }
  dl_log <- bind_rows(dl_log)
  message("Download summary: ", paste(sprintf("%s = %d", names(table(sub(" \\(.*", "", dl_log$status))),
                                              table(sub(" \\(.*", "", dl_log$status))), collapse = ", "))
}

########################################################################
## 2. Read every export in VEHSS_DIR, classify by content, de-duplicate
########################################################################

files <- c(list.files(VEHSS_DIR, pattern = "^Export.*[.]csv$", full.names = TRUE), EXTRA_FILES)
files <- files[file.exists(files)]
if (length(files) == 0) stop("No VEHSS exports found in ", VEHSS_DIR)

raw <- bind_rows(lapply(files, read_export))

file_check <- raw %>%
  group_by(source_file) %>%
  summarise(n_rows = n(),
            content = paste(unique(paste(age, sex, race, sep = " / ")), collapse = " ; "),
            n_selections = n_distinct(paste(age, sex, race)), .groups = "drop") %>%
  rowwise() %>%
  mutate(claimed = { cl <- claim_from_name(source_file)
                     if (all(is.na(cl))) NA_character_ else paste(cl["age"], cl["sex"], cl["race"], sep = " / ") }) %>%
  ungroup() %>%
  mutate(status = case_when(
    n_selections > 1   ~ "MULTIPLE selections in one file (split)",
    is.na(claimed)     ~ "name not parsed (content used)",
    claimed == content ~ "ok",
    TRUE               ~ "NAME DISAGREES WITH CONTENT (content used)"))

message("\n=== Export files: name vs content ===")
print(file_check %>% select(source_file, n_rows, claimed, content, status), n = Inf, width = Inf)
n_bad <- sum(grepl("DISAGREES", file_check$status))
if (n_bad) message(n_bad, " file(s) are labelled differently from what they contain — the content is used.")

dedup <- raw %>%
  group_by(age, sex, race) %>%
  group_modify(function(g, key) {
    fs <- unique(g$source_file)
    if (length(fs) == 1) return(g)
    wide <- g %>% select(source_file, state, prev_pct) %>%
      pivot_wider(names_from = source_file, values_from = prev_pct)
    same <- all(apply(wide[-1], 1, function(r) length(unique(r[!is.na(r)])) <= 1))
    ## prefer the downloaded file (name matches content), else the newest
    named <- fs[vapply(fs, function(f) {
      cl <- claim_from_name(f); !any(is.na(cl)) && identical(unname(cl), c(key$age, key$sex, key$race))
    }, logical(1))]
    keep <- if (length(named)) named[1] else
      g %>% group_by(source_file) %>% summarise(m = max(file_mtime), .groups = "drop") %>%
        slice_max(m, n = 1, with_ties = FALSE) %>% pull(source_file)
    if (same) message("Selection ", paste(key, collapse = " / "), " present in ", length(fs),
                      " files with identical values — keeping ", keep)
    else warning("Selection ", paste(key, collapse = " / "), " present in ", length(fs),
                 " files with DIFFERENT values (", paste(fs, collapse = ", "), ") — keeping ", keep)
    g %>% filter(source_file == keep)
  }) %>%
  ungroup()

########################################################################
## 3. Final tables
########################################################################

vehss_all <- dedup %>%
  mutate(prev = prev_pct / 100, lb = lb_pct / 100, ub = ub_pct / 100,
         pop  = cases / prev,                          # implied population of the selection
         age  = factor(age,  levels = union(AGE_LEVELS,  unique(age))),
         sex  = factor(sex,  levels = union(SEX_LEVELS,  unique(sex))),
         race = factor(race, levels = union(RACE_LEVELS, unique(race)))) %>%
  select(age, sex, race, state, state_name, year, prev, lb, ub, prev_pct, lb_pct, ub_pct,
         cases, pop, footnote, source_file) %>%
  arrange(age, sex, race, state != "US", state)

is_overall <- vehss_all$age == "40+" & vehss_all$sex == "Both" & vehss_all$race == "All"
vehss_overall <- vehss_all[is_overall, ]
vehss_combos  <- vehss_all[!is_overall, ]
if (nrow(vehss_overall) == 0)
  warning("The overall selection (40+ / Both / All races) is missing — VEHSS_comparison.R needs it.")

coverage <- expand_grid(age = AGE_LEVELS, sex = SEX_LEVELS, race = RACE_LEVELS) %>%
  left_join(vehss_all %>% group_by(age, sex, race) %>%
              summarise(n_states = sum(state != "US"), n_suppressed = sum(is.na(prev) & state != "US"),
                        has_US = any(state == "US"), source_file = first(source_file), .groups = "drop") %>%
              mutate(across(c(age, sex, race), as.character)),
            by = c("age", "sex", "race")) %>%
  mutate(present = !is.na(n_states)) %>%
  arrange(desc(present), age, sex, race)

message("\n=== Coverage: ", sum(coverage$present), " of ", nrow(coverage),
        " Age x Sex x Race selections present; suppressed state values: ",
        sum(coverage$n_suppressed, na.rm = TRUE), " ===")
miss <- coverage %>% filter(!present)
if (nrow(miss)) { message("Missing selections:"); print(miss %>% select(age, sex, race), n = Inf) }

chk <- vehss_all %>% group_by(age, sex, race) %>%
  summarise(us = sum(cases[state == "US"], na.rm = TRUE), states = sum(cases[state != "US"], na.rm = TRUE),
            rel_diff = abs(states - us) / us, .groups = "drop") %>%
  filter(is.finite(rel_diff), rel_diff > 0.01)
if (nrow(chk)) { message("Selections whose state cases do not sum to the US cases (>1%; suppressed states missing):"); print(chk, n = Inf) }

for (d in OUT_DIRS) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  write_excel_csv(vehss_overall, file.path(d, "VEHSS_AMD_overall.csv"))
  write_excel_csv(vehss_combos,  file.path(d, "VEHSS_AMD_combinations.csv"))
  write_excel_csv(coverage,      file.path(d, "VEHSS_AMD_coverage.csv"))
  message("Wrote VEHSS_AMD_overall.csv (", nrow(vehss_overall), " rows), VEHSS_AMD_combinations.csv (",
          nrow(vehss_combos), " rows), VEHSS_AMD_coverage.csv -> ", d)
}
