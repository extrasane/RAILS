#### 10_State_Prevalence_Maps.R — State-level prevalence-ratio maps
#### Run inside the AoU Researcher Workbench (the SUBGROUP workspace), AFTER:
####   Build_Disease_Indicators.R (run in the BigQuery/V9 workspace, then copy
####     its output here) -> raking_wts_w_diseases_2_code_requirement.csv
####     (person_id, state, region, w_rails, one 0/1 column per disease)
####   09_Sub_RAILS_region.R -> dt_sub_aou_region.csv (person_id, w_subrails)
#### This script merges the two by person_id (no BigQuery needed).
#### All inputs/outputs are read from / written to ../data/.
####
#### For every ICD-determined phenotype it computes each state's UNWEIGHTED
#### (AoU raw), GLOBAL-RAILS / G-RAILS (w_rails), and SUBGROUP-RAILS / S-RAILS
#### (w_subrails) prevalence, and — following the manuscript's application
#### section — the local Jensen-Shannon divergence between G-RAILS and S-RAILS
#### per state, its population-weighted total D_j, and its dispersion entropy
#### E_j (z = state). Five phenotype sets are mapped, each with the same 3
#### state-level ratio maps (G/unw, S/unw, S/G): most common (highest unweighted
#### prevalence), and the highest / lowest of D_j and of E_j — since a large D_j
#### can be misleading without knowing how concentrated (E_j) it is.
#### Graphs -> ../graph/.
####
#### AoU DISCLOSURE POLICY: any state x phenotype cell with fewer than 20 cases
#### is SUPPRESSED (shown grey / dropped), since a choropleth is an output that
#### leaves the Workbench. Contiguous-US map only: AK, HI, DC do not appear.
####
#### SEX-SPECIFIC PHENOTYPES: EXCLUDE_SEX_SPECIFIC (Settings) drops phenotypes
#### that occur in only one sex, since all prevalences use the full-population
#### denominator. Outputs of that run carry the "_nosex" suffix, so the original
#### all-phenotype outputs (EXCLUDE_SEX_SPECIFIC <- FALSE, no suffix) are kept.

library(tidyverse)
library(ggplot2)

########################################################################
## Settings
########################################################################

SUPPRESS_MIN <- 20     # AoU small-cell threshold (cases per state x phenotype)
N_SELECT     <- 12     # number of phenotypes to map (fits a facet grid)
MIN_STATES   <- 10     # a phenotype must have >= this many un-suppressed states

## rho_a (state population share, z = state). Preferred source: a by-state PUMS
## file with a person-weight column, summed to state totals. If the file is
## absent, rho_a falls back to the global-RAILS weight share (a rough proxy).
STATE_POP_FILE   <- "../data/PUMS_2022_bystate.csv"  # columns: state + a weight
STATE_POP_WEIGHT <- "PWGTP"                          # person-weight column name

## ---- Sex-specific phenotypes (reviewer comment) ----------------------------
## Every prevalence here is computed over the full adult population, not within
## the relevant sex, so a population-wide rate for a condition that occurs in
## only one sex (no BPH in females, no genital prolapse in males, ...) understates
## its condition-specific prevalence by roughly half and is not comparable to the
## other phenotypes. EXCLUDE_SEX_SPECIFIC switches the exclusion on/off; both
## versions can be kept because OUTPUT_TAG is appended to every file this script
## writes (maps, scatter, boxplots, csv tables).
##
## RULE (anatomical, not prevalence skew). A phenotype is SEX-SPECIFIC when its
## case definition can only be met by one sex, i.e. every ICD code range that
## defines it refers to (a) an organ present in only one sex - the male genital
## organs, or the female genital organs and pelvic reproductive tract - or (b) a
## state possible in only one sex - pregnancy, childbirth and the puerperium.
## Sex-skewed conditions that occur in both sexes (gout, autoimmune disease,
## osteoporosis, breast cancer, ...) are NOT sex-specific under this rule.
##
## Two implementations, chosen by SEX_SPECIFIC_RULE:
##   "icd"  : each ICD-10-CM / ICD-9-CM code range of the cause (ICD_CODE_FILE,
##            the same lookup Build_Disease_Indicators.R used) is matched, at
##            full code granularity, against the sex-specific blocks below, which
##            enumerate every category or sub-code block of ICD-10-CM / ICD-9-CM
##            that is restricted to one sex (genital organs, pregnancy). A cause is
##            sex-specific iff ALL of its DIAGNOSTIC ranges match AND all matches
##            are the SAME sex. Supplementary codes (ICD-10 Z, ICD-9 V/E: screening,
##            personal / family history, status, external cause) describe an
##            encounter rather than a disease and are ignored in the classification
##            when IGNORE_SUPPLEMENTARY_CODES is TRUE - otherwise e.g. prostate
##            cancer, whose GBD code list carries Z12.5 / Z85.46 / Z80.42, would be
##            called "mixed". A cause whose diagnostic ranges are only partly
##            sex-specific, or that combines male and female blocks (an aggregate
##            like "Urogenital congenital anomalies" or "Urinary diseases and male
##            infertility"), is "mixed" and is KEPT. Every cause's class and its
##            per-range sexes are written to a csv so mixed causes can be reviewed
##            by hand; hand decisions go in SEX_SPECIFIC_MANUAL_EXCLUDE / _KEEP.
##   "name" : case-insensitive substrings of the cause name; used automatically
##            when ICD_CODE_FILE is not available in this workspace.
## Breast cancer occurs in both sexes (~1% male) so the strict rule keeps it; set
## EXCLUDE_BREAST_CANCER <- TRUE to drop it as well under either rule.
EXCLUDE_SEX_SPECIFIC       <- TRUE       # TRUE: drop them | FALSE: original full run
SEX_SPECIFIC_RULE          <- "icd"      # "icd" or "name"
EXCLUDE_BREAST_CANCER      <- FALSE
IGNORE_SUPPLEMENTARY_CODES <- TRUE       # drop ICD-10 Z / ICD-9 V, E ranges before classifying
ICD_CODE_FILE              <- "../data/all_icd_codes.csv"   # ICD9CM / ICD10CM / Cause
OUTPUT_TAG                 <- if (EXCLUDE_SEX_SPECIFIC) "_nosex" else ""

## Hand overrides applied AFTER the rule (exact cause names). Use them for
## aggregates the rule calls "mixed" but that are one-sex in an adult cohort, e.g.
## "Maternal and neonatal disorders" (neonatal codes cannot occur in AoU adults,
## so every case is maternal), or to keep a cause the rule excluded.
SEX_SPECIFIC_MANUAL_EXCLUDE <- c(
  ## aggregate incl. neonatal P codes, which cannot occur in AoU adults
  "Maternal and neonatal disorders",
  ## every range is F except the GBD range "F53-F54", which runs past the
  ## puerperal category F53 into F54 (psychological factors in diseases
  ## classified elsewhere) - a mapping artefact, not a both-sex condition
  "Maternal disorders", "Other maternal disorders",
  ## the GBD definition also carries breast disorders (N61-N64 / 611), urogenital
  ## candidiasis (B37.4 / 112.2), urinary symptoms (R30-R39 / 788) and decreased
  ## libido (799.81), all of which men can have, so the strict rule says "mixed";
  ## the cause is female-only in intent and name, and a man carrying only those
  ## codes would be a phenotyping artefact rather than a gynecological case
  "Gynecological diseases", "Other gynecological diseases",
  ## GBD aggregate that lumps a both-sex group (urinary diseases) with a
  ## male-only one (male infertility); excluded on clinical advice because the
  ## lumped definition may be sex-specific in part
  "Urinary diseases and male infertility"
)
## Reviewed and KEPT although "mixed": "Benign and in situ intestinal neoplasms"
## (its GBD list contains 236.0 = uncertain-behaviour neoplasm of uterus, an
## apparent mapping error), "Other chromosomal abnormalities" (contains Q97 with
## non-sex-specific Q9x codes) and every aggregate cause spanning both sexes.
SEX_SPECIFIC_MANUAL_KEEP    <- character(0)

## ---- Sex-specific ICD blocks -------------------------------------------------
## Codes are compared with dots removed, as prefix intervals: a range [f, l] lies
## in a block [s, e] iff f >= s and l is not beyond e (so "D07" is NOT inside
## "D070-D073", but "D075" IS inside "D074-D076"; "O9A" > "O99"). Sub-code
## entries (4-5 characters) cover the sex-specific parts of mixed categories.
ICD10_SEX_BLOCKS <- tibble::tribble(
  ~start,  ~end,    ~sex, ~what,
  ## ---- male ----
  "N40",   "N53",   "M",  "diseases of male genital organs (BPH N40, male infertility N46, erectile/sexual dysfunction N52-N53)",
  "C60",   "C63",   "M",  "malignant neoplasms of male genital organs (penis, prostate, testis, other)",
  "D074",  "D076",  "M",  "carcinoma in situ of penis / prostate / other male genital organs",
  "D29",   "D29",   "M",  "benign neoplasm of male genital organs",
  "D40",   "D40",   "M",  "neoplasm of uncertain behaviour of male genital organs",
  "E29",   "E29",   "M",  "testicular dysfunction",
  "L291",  "L291",  "M",  "pruritus scroti",
  "Q53",   "Q55",   "M",  "congenital malformations of male genital organs (undescended testis, hypospadias, other)",
  "Q98",   "Q98",   "M",  "other sex chromosome abnormalities, male phenotype (Klinefelter etc.)",
  "R86",   "R86",   "M",  "abnormal findings in specimens from male genital organs",
  "Z125",  "Z125",  "M",  "screening for malignant neoplasm of prostate",
  "Z1271", "Z1271", "M",  "screening for malignant neoplasm of testis",
  "Z8042", "Z8043", "M",  "family history of malignant neoplasm of prostate / testis",
  "Z8545", "Z8548", "M",  "personal history of malignant neoplasm of male genital organs",
  "Z8743", "Z8743", "M",  "personal history of diseases of male genital organs",
  ## ---- female ----
  "N70",   "N98",   "F",  "diseases of female pelvic organs / genital tract (PID N70-N77, endometriosis N80, prolapse N81, menstrual/menopausal N91-N95, female infertility N97)",
  "N993",  "N993",  "F",  "prolapse of vaginal vault after hysterectomy",
  "C51",   "C58",   "F",  "malignant neoplasms of female genital organs (vulva, vagina, cervix, uterus, ovary, other, placenta)",
  "D06",   "D06",   "F",  "carcinoma in situ of cervix uteri",
  "D070",  "D073",  "F",  "carcinoma in situ of endometrium / vulva / vagina / other female genital organs",
  "D25",   "D28",   "F",  "benign neoplasms of female genital organs (leiomyoma, uterus, ovary, other)",
  "D39",   "D39",   "F",  "neoplasm of uncertain behaviour of female genital organs",
  "E28",   "E28",   "F",  "ovarian dysfunction (incl. polycystic ovarian syndrome E28.2)",
  "B373",  "B373",  "F",  "candidiasis of vulva and vagina",
  "L292",  "L292",  "F",  "pruritus vulvae",
  "Q50",   "Q52",   "F",  "congenital malformations of female genital organs",
  "Q96",   "Q97",   "F",  "Turner syndrome; other sex chromosome abnormalities, female phenotype",
  "R87",   "R87",   "F",  "abnormal findings in specimens from female genital organs",
  "A34",   "A34",   "F",  "obstetrical tetanus",
  "F53",   "F53",   "F",  "mental and behavioural disorders associated with the puerperium",
  "M830",  "M830",  "F",  "puerperal osteomalacia",
  "O00",   "O9A",   "F",  "pregnancy, childbirth and the puerperium (whole chapter)",
  "Z124",  "Z124",  "F",  "screening for malignant neoplasm of cervix",
  "Z1272", "Z1273", "F",  "screening for malignant neoplasm of vagina / ovary",
  "Z32",   "Z37",   "F",  "encounters: pregnancy test, pregnant state, supervision of pregnancy, antenatal screening, outcome of delivery",
  "Z39",   "Z39",   "F",  "encounter for postpartum care",
  "Z3A",   "Z3A",   "F",  "weeks of gestation",
  "Z640",  "Z641",  "F",  "problems related to unwanted pregnancy / multiparity",
  "Z8041", "Z8041", "F",  "family history of malignant neoplasm of ovary",
  "Z8540", "Z8544", "F",  "personal history of malignant neoplasm of female genital organs",
  "Z8741", "Z8742", "F",  "personal history of cervical dysplasia / other female genital diseases",
  "Z875",  "Z875",  "F",  "personal history of complications of pregnancy, childbirth and the puerperium",
  "Z9071", "Z9072", "F",  "acquired absence of cervix and uterus / ovaries",
  "Z975",  "Z975",  "F",  "presence of (intrauterine) contraceptive device",
  "Z98891","Z98891","F",  "history of uterine scar from previous surgery"
)
## ICD-9-CM equivalents
ICD9_SEX_BLOCKS <- tibble::tribble(
  ~start,  ~end,    ~sex, ~what,
  ## ---- male ----
  "600",   "608",   "M",  "diseases of male genital organs (BPH 600, male infertility 606, impotence 607.84 within 607)",
  "185",   "187",   "M",  "malignant neoplasms of male genital organs",
  "2334",  "2336",  "M",  "carcinoma in situ of prostate / penis / other male genital organs",
  "222",   "222",   "M",  "benign neoplasm of male genital organs",
  "2364",  "2366",  "M",  "neoplasm of uncertain behaviour of testis / prostate / other male genital organs",
  "257",   "257",   "M",  "testicular dysfunction",
  "7525",  "7526",  "M",  "undescended testis; hypospadias and epispadias",
  "7587",  "7587",  "M",  "Klinefelter's syndrome",
  "V1045", "V1049", "M",  "personal history of malignant neoplasm of male genital organs",
  "V1642", "V1643", "M",  "family history of malignant neoplasm of prostate / testis",
  "V7644", "V7645", "M",  "screening for malignant neoplasm of prostate / testis",
  ## ---- female ----
  "614",   "629",   "F",  "disorders of female genital tract (PID 614-616, endometriosis 617, prolapse 618, menstrual/menopausal 625-627, female infertility 628)",
  "179",   "184",   "F",  "malignant neoplasms of female genital organs",
  "2331",  "2333",  "F",  "carcinoma in situ of cervix / other and unspecified uterus / other female genital organs",
  "218",   "221",   "F",  "benign neoplasms of female genital organs",
  "2360",  "2363",  "F",  "neoplasm of uncertain behaviour of uterus / placenta / ovary / other female genital organs",
  "256",   "256",   "F",  "ovarian dysfunction",
  "1121",  "1121",  "F",  "candidiasis of vulva and vagina",
  "7520",  "7524",  "F",  "congenital anomalies of ovaries, fallopian tubes, uterus, cervix, vagina",
  "7586",  "7586",  "F",  "gonadal dysgenesis (Turner's syndrome)",
  "7950",  "7951",  "F",  "abnormal Papanicolaou smear of cervix / vagina",
  "630",   "679",   "F",  "complications of pregnancy, childbirth and the puerperium (whole chapter)",
  "V132",  "V132",  "F",  "personal history of other genital system and obstetric disorders",
  "V1040", "V1044", "F",  "personal history of malignant neoplasm of female genital organs",
  "V1641", "V1641", "F",  "family history of malignant neoplasm of ovary",
  "V22",   "V24",   "F",  "supervision of pregnancy / postpartum care",
  "V27",   "V28",   "F",  "outcome of delivery / antenatal screening",
  "V615",  "V617",  "F",  "multiparity; illegitimacy; unwanted pregnancy",
  "V723",  "V724",  "F",  "gynecological examination; pregnancy examination or test",
  "V762",  "V762",  "F",  "screening for malignant neoplasm of cervix",
  "V7646", "V7647", "F",  "screening for malignant neoplasm of ovary / vagina"
)
if (EXCLUDE_BREAST_CANCER) {
  ICD10_SEX_BLOCKS <- dplyr::add_row(ICD10_SEX_BLOCKS, start = "C50", end = "C50",
                                     sex = "F", what = "malignant neoplasm of breast (by convention)")
  ICD9_SEX_BLOCKS  <- dplyr::add_row(ICD9_SEX_BLOCKS,  start = "174", end = "175",
                                     sex = "F", what = "malignant neoplasm of breast (by convention)")
}

## Name patterns for SEX_SPECIFIC_RULE == "name" (GBD-style cause names). Kept
## deliberately specific ("cervical cancer", not "cervical", which would also hit
## cervical-spine conditions).
SEX_SPECIFIC_PATTERNS <- c(
  ## male
  "prostat", "testic", "erectile", "penile", "hypospadias", "undescended test",
  "male infertility",                       # also matches "female infertility"
  ## female
  "genital prolapse", "uterine prolapse", "pelvic organ prolapse",
  "endometriosis", "polycystic ovar", "ovarian", "ovary", "uterine", "uterus",
  "vagin", "vulv", "cervix", "pelvic inflammatory", "pelvic organ", "menstrua",
  "menopaus", "gynecolog", "cervical cancer",
  ## pregnancy   ("genital" alone is NOT used: it would match "congenital")
  "maternal", "pregnan", "ectopic", "puerper", "postpartum", "obstetric",
  "gestation", "abortion", "miscarriage",
  if (EXCLUDE_BREAST_CANCER) "breast cancer"
)

## Split "A00-A09, B10" style code lists into (first, last) category ranges
split_ranges <- function(code_string) {
  if (is.na(code_string) || !nzchar(trimws(code_string))) return(NULL)
  parts <- trimws(strsplit(code_string, ",")[[1]])
  parts <- parts[nzchar(parts)]
  fl <- strsplit(parts, "-")
  tibble::tibble(first = vapply(fl, `[`, "", 1),
                 last  = vapply(fl, function(x) if (length(x) > 1) x[2] else x[1], ""))
}
## Normalise a code for comparison: upper-case, no dots / spaces ("D07.5" -> "D075")
norm_code <- function(x) toupper(gsub("[[:space:].]", "", x))

## Supplementary (non-diagnostic) ranges: ICD-10 Z chapter; ICD-9 V and E codes
is_supplementary <- function(first, icd) {
  f <- norm_code(first)
  if (icd == 10) startsWith(f, "Z") else grepl("^[VE]", f)
}

## Sex ("M"/"F") of the block containing each range, NA if it lies in no block.
## Prefix-interval containment: f >= s, and l (as a prefix interval [l, l~]) does
## not extend past e ('~' sorts after every digit and letter IN BYTE ORDER, so the
## comparisons are done under the C collation - the session locale may sort
## punctuation first and would silently break the sub-code blocks).
range_sex <- function(rng, blocks) {
  if (is.null(rng) || nrow(rng) == 0) return(character(0))
  old_collate <- Sys.getlocale("LC_COLLATE")
  invisible(Sys.setlocale("LC_COLLATE", "C"))
  on.exit(invisible(Sys.setlocale("LC_COLLATE", old_collate)), add = TRUE)
  vapply(seq_len(nrow(rng)), function(i) {
    f <- norm_code(rng$first[i]); l <- norm_code(rng$last[i])
    hit <- blocks$sex[f >= blocks$start & paste0(l, "~") <= paste0(blocks$end, "~")]
    if (length(hit)) hit[1] else NA_character_
  }, character(1))
}

## Per-range sexes of one cause (diagnostic ranges only when
## IGNORE_SUPPLEMENTARY_CODES), as a named character vector for the csv
cause_range_sexes <- function(icd10, icd9) {
  r10 <- split_ranges(icd10); r9 <- split_ranges(icd9)
  if (IGNORE_SUPPLEMENTARY_CODES) {
    if (!is.null(r10)) r10 <- r10[!vapply(r10$first, is_supplementary, logical(1), icd = 10), , drop = FALSE]
    if (!is.null(r9))  r9  <- r9[!vapply(r9$first,  is_supplementary, logical(1), icd = 9),  , drop = FALSE]
  }
  s <- c(range_sex(r10, ICD10_SEX_BLOCKS), range_sex(r9, ICD9_SEX_BLOCKS))
  lab <- c(if (!is.null(r10) && nrow(r10)) paste0("10:", r10$first, ifelse(r10$last != r10$first, paste0("-", r10$last), "")),
           if (!is.null(r9)  && nrow(r9))  paste0("9:",  r9$first,  ifelse(r9$last  != r9$first,  paste0("-", r9$last),  "")))
  setNames(s, lab)
}

## Classify one cause: "M" / "F" (all diagnostic ranges in one sex's blocks),
## "mixed" (some ranges sex-specific, or male and female blocks combined), "none".
classify_cause <- function(icd10, icd9) {
  s <- cause_range_sexes(icd10, icd9)
  if (length(s) == 0 || all(is.na(s))) return("none")
  if (any(is.na(s)))                    return("mixed")
  if (length(unique(s)) > 1)            return("mixed")
  unique(unname(s))
}

## "icd" rule. Returns the sex-specific cause names, and writes the full
## classification (class, per-range sexes, code ranges) so mixed causes can be
## reviewed. Returns NULL if the lookup file is missing (caller falls back to "name").
sex_specific_by_icd <- function(disease_names) {
  if (!file.exists(ICD_CODE_FILE)) {
    warning("ICD_CODE_FILE not found (", ICD_CODE_FILE, ") - falling back to the name rule.")
    return(NULL)
  }
  codes <- readr::read_csv(ICD_CODE_FILE, show_col_types = FALSE) %>%
    dplyr::filter(!(is.na(ICD10CM) & is.na(ICD9CM)))
  cls <- codes %>%
    dplyr::mutate(sex_class   = mapply(classify_cause, ICD10CM, ICD9CM),
                  range_sexes = mapply(function(a, b) {
                    s <- cause_range_sexes(a, b)
                    paste0(names(s), "=", ifelse(is.na(s), "-", s), collapse = "; ")
                  }, ICD10CM, ICD9CM),
                  analysed    = Cause %in% disease_names) %>%
    dplyr::select(Cause, sex_class, analysed, range_sexes, ICD9CM, ICD10CM) %>%
    dplyr::arrange(sex_class, Cause)
  readr::write_excel_csv(cls, paste0("../data/sex_specific_icd_classification", OUTPUT_TAG, ".csv"))
  mixed <- cls$Cause[cls$sex_class == "mixed" & cls$analysed]
  if (length(mixed))
    message("Causes with PARTLY sex-specific codes (kept; review by hand): ",
            paste(mixed, collapse = "; "))
  not_in_lookup <- setdiff(disease_names, codes$Cause)
  if (length(not_in_lookup))
    warning("Analysed phenotypes not found in ICD_CODE_FILE (kept, unclassified): ",
            paste(not_in_lookup, collapse = "; "))
  intersect(cls$Cause[cls$sex_class %in% c("M", "F")], disease_names)
}

## "name" rule
sex_specific_by_name <- function(disease_names) {
  pats <- SEX_SPECIFIC_PATTERNS[!vapply(SEX_SPECIFIC_PATTERNS, is.null, logical(1))]
  disease_names[Reduce(`|`, lapply(pats, function(p)
    grepl(p, disease_names, ignore.case = TRUE)))]
}

########################################################################
## Load the disease matrix (person_id, state, region, w_rails, <diseases>)
## and MERGE the regional sub-RAILS weight (w_subrails) by person_id.
## Disease columns are everything that is not an id / state / weight column.
########################################################################

dz   <- read_csv("../data/raking_wts_w_diseases_2_code_requirement.csv",
                 show_col_types = FALSE)
subw <- read_csv("../data/dt_sub_aou_region.csv", show_col_types = FALSE) %>%
  select(person_id, w_subrails)

df0 <- dz %>% left_join(subw, by = "person_id")

id_wt_cols   <- c("person_id", "state", "region", "w_rails", "w_subrails")
disease_cols <- setdiff(names(df0), id_wt_cols)

## Sex-specific phenotypes: excluded only when EXCLUDE_SEX_SPECIFIC is TRUE
## (settings above). With FALSE this block is a no-op and the run reproduces the
## original, all-phenotype outputs (OUTPUT_TAG = "").
sex_specific <- character(0)
if (EXCLUDE_SEX_SPECIFIC) {
  rule_used <- SEX_SPECIFIC_RULE
  sex_specific <- switch(SEX_SPECIFIC_RULE,
                         icd  = sex_specific_by_icd(disease_cols),
                         name = sex_specific_by_name(disease_cols),
                         stop("SEX_SPECIFIC_RULE must be \"icd\" or \"name\""))
  if (is.null(sex_specific)) {                 # lookup file missing -> name rule
    rule_used    <- "name"
    sex_specific <- sex_specific_by_name(disease_cols)
  }
  message("Sex-specific phenotypes excluded by the '", rule_used, "' rule (",
          length(sex_specific), "): ", paste(sex_specific, collapse = "; "))
  ## hand overrides (exact names; only those actually present are applied)
  man_ex <- intersect(SEX_SPECIFIC_MANUAL_EXCLUDE, disease_cols)
  man_kp <- intersect(SEX_SPECIFIC_MANUAL_KEEP, sex_specific)
  if (length(man_ex)) message("Manually excluded in addition: ", paste(man_ex, collapse = "; "))
  if (length(man_kp)) message("Manually kept despite the rule: ", paste(man_kp, collapse = "; "))
  excl_tab <- dplyr::bind_rows(
    tibble::tibble(Disease = setdiff(sex_specific, man_kp), rule = rule_used),
    tibble::tibble(Disease = setdiff(man_ex, sex_specific), rule = "manual"))
  sex_specific <- excl_tab$Disease
  disease_cols <- setdiff(disease_cols, sex_specific)
  write_excel_csv(excl_tab, paste0("../data/sex_specific_phenotypes_excluded", OUTPUT_TAG, ".csv"))
}
message("Phenotypes analysed: ", length(disease_cols),
        if (EXCLUDE_SEX_SPECIFIC) " (after the sex-specific exclusion)" else "")

df <- df0 %>% filter(!is.na(state), !is.na(w_rails), !is.na(w_subrails))
message("Participants with state + both weights: ", nrow(df),
        " | ", length(disease_cols), " diseases")

########################################################################
## Three state-level prevalences per phenotype: UNWEIGHTED (AoU raw),
## GLOBAL RAILS (w_rails), and REGIONAL sub-RAILS (w_subrails). The maps show
## the pairwise RAW ratios between them (0 to infinity; 1 = equal).
##
## Computed with SPARSE MATRIX PRODUCTS rather than a pivot_longer melt:
## an N-participant x ~300-disease table would explode into a ~100M-row frame
## and kill the kernel. Here each quantity is a single crossprod().
########################################################################

library(Matrix)

D    <- as(data.matrix(df[disease_cols]), "dgCMatrix")  # persons x diseases (0/1)
sfac <- factor(df$state)
S    <- sparse.model.matrix(~ 0 + sfac)                 # persons x states
states <- levels(sfac)
w_g  <- df$w_rails
w_r  <- df$w_subrails

## National prevalence per disease under each scheme (unweighted / G-RAILS / S-RAILS)
P_unw   <- setNames(as.numeric(Matrix::colSums(D)) / nrow(df), disease_cols)
P_g_nat <- setNames(as.numeric(crossprod(D, w_g)) / sum(w_g), disease_cols)
P_r_nat <- setNames(as.numeric(crossprod(D, w_r)) / sum(w_r), disease_cols)

## Per state x disease: weighted case sums, raw case counts, weight/size totals
num_g   <- as.matrix(crossprod(S, Diagonal(x = w_g) %*% D))   # states x diseases
num_r   <- as.matrix(crossprod(S, Diagonal(x = w_r) %*% D))
ncase   <- as.matrix(crossprod(S, D))
tot_g   <- as.numeric(crossprod(S, w_g))                      # per state weight totals
tot_r   <- as.numeric(crossprod(S, w_r))
n_state <- as.numeric(Matrix::colSums(S))                     # per state participant count
dimnames(num_g) <- dimnames(num_r) <- dimnames(ncase) <- list(states, disease_cols)
rm(D, S); gc()

########################################################################
## Application metrics (manuscript application section), with z = STATE.
##   rho_a          : state population share N^a/N (from STATE_POP_FILE if
##                    present, else the global-RAILS weight share as a proxy)
##   d^a_{j,GS}      : local Jensen-Shannon divergence between the G-RAILS and
##                    S-RAILS Bernoulli prevalences in state a
##   D_j = sum_a rho_a d^a_{j,GS}          : population-weighted G-vs-S divergence
##   E_j = -sum_a q_a log q_a / log L      : how evenly D_j is spread across states
##         (q_a = rho_a d^a / D_j; 0 = one state dominates, 1 = uniform)
########################################################################

rho_vec <- local({
  if (file.exists(STATE_POP_FILE)) {
    sp <- readr::read_csv(STATE_POP_FILE, show_col_types = FALSE)
    if (all(c("state", STATE_POP_WEIGHT) %in% names(sp))) {
      pop <- sp %>% group_by(state) %>%
        summarise(pop = sum(.data[[STATE_POP_WEIGHT]], na.rm = TRUE), .groups = "drop")
      v <- setNames(pop$pop, pop$state)[states]
      v[is.na(v)] <- 0
      message("rho_a from ", STATE_POP_FILE, " (weight = ", STATE_POP_WEIGHT, ")")
      return(setNames(v / sum(v), states))
    }
    warning("STATE_POP_FILE lacks 'state' / '", STATE_POP_WEIGHT,
            "' — using global-RAILS weight share for rho_a.")
  } else {
    message("STATE_POP_FILE not found — rho_a from global-RAILS weight share (proxy).")
  }
  setNames(tot_g / sum(w_g), states)   # fallback proxy
})

H_bern  <- function(p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); -p * log(p) - (1 - p) * log(1 - p) }
js_bern <- function(pg, ps) H_bern((pg + ps) / 2) - 0.5 * H_bern(pg) - 0.5 * H_bern(ps)

weighted_js  <- function(rho, d) { ok <- is.finite(d) & is.finite(rho); sum(rho[ok] * d[ok]) }
disc_entropy <- function(rho, d) {
  ok  <- is.finite(d) & is.finite(rho)
  c_a <- rho[ok] * d[ok]
  Dj  <- sum(c_a)
  if (Dj <= 0)           return(NA_real_)   # G-RAILS == S-RAILS everywhere
  if (length(c_a) <= 1)  return(0)
  q <- c_a / Dj; q <- q[q > 0]
  -sum(q * log(q)) / log(length(c_a))
}

## Melt the small (states x diseases) matrices to long form
melt_mat <- function(m, value) {
  as.data.frame(m) %>%
    tibble::rownames_to_column("state") %>%
    pivot_longer(-state, names_to = "Disease", values_to = value)
}

## Three state-level prevalences per phenotype:
##   p_unw = unweighted (mu_N), p_g = G-RAILS (mu_G), p_r = S-RAILS (mu_S)
state_tab <- melt_mat(ncase, "n_cases") %>%
  left_join(melt_mat(sweep(ncase, 1, n_state, "/"), "p_unw"), by = c("state", "Disease")) %>%
  left_join(melt_mat(sweep(num_g, 1, tot_g,   "/"), "p_g"),   by = c("state", "Disease")) %>%
  left_join(melt_mat(sweep(num_r, 1, tot_r,   "/"), "p_r"),   by = c("state", "Disease")) %>%
  mutate(
    suppress = n_cases < SUPPRESS_MIN,
    rho      = rho_vec[state],
    ## Prevalence ratios (0 to infinity; 1 = equal)
    r_g_unw = ifelse(suppress | p_unw <= 0, NA_real_, p_g / p_unw),  # G-RAILS vs unweighted
    r_s_unw = ifelse(suppress | p_unw <= 0, NA_real_, p_r / p_unw),  # S-RAILS vs unweighted
    r_s_g   = ifelse(suppress | p_g   <= 0, NA_real_, p_r / p_g),    # S-RAILS vs G-RAILS (R_{S/G})
    ## Local Jensen-Shannon divergence between G-RAILS and S-RAILS (d^a_{j,GS})
    d_gs    = ifelse(suppress, NA_real_, js_bern(p_g, p_r))
  )

message("Suppressed state x phenotype cells (<", SUPPRESS_MIN, " cases): ",
        sum(state_tab$suppress))

########################################################################
## Phenotype selection — each set mapped with the same 3 state-level ratios.
## A large D_j can be misleading without its dispersion E_j, so map the extremes
## of BOTH: most common, highest/lowest D_j, and highest/lowest E_j.
########################################################################

pheno <- state_tab %>%
  group_by(Disease) %>%
  summarise(
    n_states         = sum(!suppress & is.finite(d_gs)),
    D_js             = weighted_js(rho, d_gs),   # population-weighted G-vs-S divergence
    E_js             = disc_entropy(rho, d_gs),  # how evenly D_j is spread across states
    sum_abs_logratio = sum(abs(log(r_s_g)), na.rm = TRUE),  # total |log(S/G)| over states
    .groups  = "drop"
  ) %>%
  mutate(P_unw   = P_unw[Disease],
         P_g_nat = P_g_nat[Disease],
         P_r_nat = P_r_nat[Disease]) %>%
  filter(n_states >= MIN_STATES)
message("Phenotypes with >= ", MIN_STATES, " un-suppressed states (the manuscript's ",
        "phenotype count): ", nrow(pheno))

sel_hi <- function(col) pheno %>% filter(is.finite(.data[[col]])) %>%
  slice_max(.data[[col]], n = N_SELECT, with_ties = FALSE) %>% pull(Disease)
sel_lo <- function(col) pheno %>% filter(is.finite(.data[[col]])) %>%
  slice_min(.data[[col]], n = N_SELECT, with_ties = FALSE) %>% pull(Disease)

setA     <- pheno %>% slice_max(P_unw, n = N_SELECT, with_ties = FALSE) %>% pull(Disease)
setHighD <- sel_hi("D_js")     # largest magnitude of G-vs-S divergence
setLowD  <- sel_lo("D_js")     # smallest
setHighE <- sel_hi("E_js")     # divergence spread most evenly across states
setLowE  <- sel_lo("E_js")     # divergence most concentrated in a few states

message("Set A (most common): ", paste(setA,     collapse = "; "))
message("Highest D:           ", paste(setHighD, collapse = "; "))
message("Lowest D:            ", paste(setLowD,  collapse = "; "))
message("Highest E:           ", paste(setHighE, collapse = "; "))
message("Lowest E:            ", paste(setLowE,  collapse = "; "))

## Ranked summary: national prevalence (unweighted / G-RAILS / S-RAILS) + D + E,
## printed for the top and bottom phenotypes by D and, separately, by E.
rank_tab <- pheno %>%
  transmute(Disease,
            prev_unw      = P_unw,
            prev_global   = P_g_nat,
            prev_subgroup = P_r_nat,
            D = D_js, E = E_js)

print_rank <- function(tab, col) {
  ord <- tab %>% arrange(desc(.data[[col]]))
  message("\n--- Top ", N_SELECT, " phenotypes by ", col, " ---")
  print(head(ord, N_SELECT), n = N_SELECT)
  message("\n--- Bottom ", N_SELECT, " phenotypes by ", col, " ---")
  print(tail(ord, N_SELECT), n = N_SELECT)
}

print_rank(rank_tab, "D")
print_rank(rank_tab, "E")

########################################################################
## Representative D-E regime examples (manuscript Table "DEexamples").
## The interpretation table reads the pair (D_j, E_j) in four regimes:
##   high D, low E  : G- and S-RAILS differ, driven by one / a few states
##   high D, high E : they differ everywhere by comparable amounts
##   low D,  low E  : they agree overall; the small difference sits in one state
##   low D,  high E : they agree everywhere
## "High"/"low" are cut at quantiles of D_j and E_j over the analysed phenotypes
## (REGIME_*). Within each regime the N_EXAMPLES phenotypes closest to that corner
## (by percentile rank) are listed together with the LEVERAGE SUBGROUP: the state
## carrying the largest share q_a = rho_a d^a_j / D_j of the divergence. It is
## reported as the leverage state when that share is at least LEVERAGE_SHARE_MIN
## (a majority of D_j), otherwise the discrepancy is reported as "spread". The
## share is the only criterion; no condition is placed on the state's own ratio.
## Output: ../data/DE_regime_examples<TAG>.csv, ../graph/plot_D_vs_E_regimes<TAG>.png
########################################################################

## Defaults: "high D" = top decile, "low D" = bottom quartile (D is heavily
## right-skewed); "high"/"low" E = above/below the median, so every high-D
## phenotype falls in one of the two high-D regimes. Adjust to taste.
REGIME_D_HI     <- 0.90   # D_j at or above this quantile  -> "high D"
REGIME_D_LO     <- 0.25   # D_j at or below this quantile  -> "low D"
REGIME_E_HI     <- 0.50   # E_j at or above this quantile  -> "high E"
REGIME_E_LO     <- 0.50   # E_j below this quantile        -> "low E"
N_EXAMPLES      <- 3
LEVERAGE_SHARE_MIN <- 0.50   # leverage state = the state with the largest share q_a of
                             # D_j, reported when that share is a MAJORITY; otherwise
                             # "spread". (E is the evenness of ALL shares, so one state
                             # at ~36% with the rest spread evenly still gives a high E.)
MAP_RATIO_CAP      <- 2.5    # regime maps: symmetric log colour scale, ratios squished to
                             # [1/cap, cap] so one extreme state does not wash out the rest

REGIME_READING <- c(
  "High D, low E"  = "Global and subgroup weights disagree, and most of the divergence sits in the leverage state; explain that state disease by disease (place of diagnosis vs residence, site specialization) and inspect or set aside it before reporting state-level estimates.",
  "High D, high E" = "The calibration choice shifts prevalence in many states by comparable amounts; subgroup weights are warranted if regional selection is believed to differ.",
  "Low D, low E"   = "The two strategies agree overall; global weights suffice, with the leverage state noted if it is the analytic target.",
  "Low D, high E"  = "Global weights suffice everywhere; subgroup calibration is unnecessary."
)

## State carrying the largest share of each phenotype's divergence
top_state <- state_tab %>%
  filter(!suppress, is.finite(d_gs), is.finite(rho)) %>%
  group_by(Disease) %>%
  mutate(q = rho * d_gs / sum(rho * d_gs)) %>%
  slice_max(q, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Disease, top_state = state, top_share = q, top_state_ratio_SG = r_s_g)

regime_cuts <- with(dplyr::filter(pheno, is.finite(D_js), is.finite(E_js)), c(
  D_hi = unname(quantile(D_js, REGIME_D_HI)), D_lo = unname(quantile(D_js, REGIME_D_LO)),
  E_hi = unname(quantile(E_js, REGIME_E_HI)), E_lo = unname(quantile(E_js, REGIME_E_LO))))
message("Regime cut-points: ", paste(names(regime_cuts), signif(regime_cuts, 3),
                                     sep = " = ", collapse = "; "))

regime_all <- pheno %>%
  filter(is.finite(D_js), is.finite(E_js)) %>%
  left_join(top_state, by = "Disease") %>%
  mutate(
    pD = percent_rank(D_js), pE = percent_rank(E_js),
    regime = case_when(
      D_js >= regime_cuts["D_hi"] & E_js <  regime_cuts["E_lo"] ~ "High D, low E",
      D_js >= regime_cuts["D_hi"] & E_js >= regime_cuts["E_hi"] ~ "High D, high E",
      D_js <= regime_cuts["D_lo"] & E_js <  regime_cuts["E_lo"] ~ "Low D, low E",
      D_js <= regime_cuts["D_lo"] & E_js >= regime_cuts["E_hi"] ~ "Low D, high E",
      TRUE ~ NA_character_),
    ## distance towards the regime's corner, larger = more typical of the regime
    corner_score = case_when(
      regime == "High D, low E"  ~ pD + (1 - pE),
      regime == "High D, high E" ~ pD + pE,
      regime == "Low D, low E"   ~ (1 - pD) + (1 - pE),
      regime == "Low D, high E"  ~ (1 - pD) + pE),
    leverage_state = ifelse(!is.na(top_share) & top_share >= LEVERAGE_SHARE_MIN,
                            top_state, "spread"))

regime_examples <- regime_all %>%
  filter(!is.na(regime)) %>%
  group_by(regime) %>%
  slice_max(corner_score, n = N_EXAMPLES, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(regime = factor(regime, levels = names(REGIME_READING)),
         practical_reading = REGIME_READING[as.character(regime)]) %>%
  arrange(regime, desc(corner_score)) %>%
  select(regime, Disease, prev_unw = P_unw, prev_global = P_g_nat, prev_subgroup = P_r_nat,
         D = D_js, E = E_js, leverage_state, top_state, top_share, top_state_ratio_SG,
         practical_reading)

message("\n--- Representative D-E regime examples ---")
print(regime_examples %>% select(-practical_reading) %>%
        mutate(D = signif(D, 3), E = round(E, 3), top_share = round(top_share, 2)),
      n = Inf, width = Inf)
write_excel_csv(regime_examples, paste0("../data/DE_regime_examples", OUTPUT_TAG, ".csv"))

## D-E scatter with the regime cut-points and the selected examples labelled
regime_cols <- c("High D, low E" = "#B2182B", "High D, high E" = "#E08214",
                 "Low D, low E" = "#2166AC", "Low D, high E" = "#1B7837")
p_regime <- ggplot(regime_all, aes(D_js, E_js)) +
  geom_point(colour = "grey80", size = 1.6, alpha = 0.8) +
  geom_vline(xintercept = regime_cuts[c("D_lo", "D_hi")], linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = regime_cuts[c("E_lo", "E_hi")], linetype = "dashed", colour = "grey50") +
  geom_point(data = regime_examples, aes(D, E, colour = regime), size = 3) +
  scale_colour_manual(values = regime_cols, name = NULL) +
  labs(x = expression(D[j]), y = expression(E[j]),
       title = "D-E regimes: representative phenotypes",
       subtitle = sprintf("dashed lines: D quantiles %.2f / %.2f, E quantiles %.2f / %.2f",
                          REGIME_D_LO, REGIME_D_HI, REGIME_E_LO, REGIME_E_HI)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_regime <- p_regime + ggrepel::geom_text_repel(
    data = regime_examples,
    aes(D, E, label = ifelse(leverage_state == "spread", Disease,
                             paste0(Disease, " (", leverage_state, ")")),
        colour = regime),
    size = 3, box.padding = 0.6, point.padding = 0.3, min.segment.length = 0,
    segment.size = 0.3, segment.color = "grey60", max.overlaps = Inf, seed = 1,
    show.legend = FALSE)
}
print(p_regime)
ggsave(paste0("../graph/plot_D_vs_E_regimes", OUTPUT_TAG, ".png"), p_regime,
       width = 9.5, height = 7, dpi = 300)

########################################################################
## US state polygons (contiguous US; from the `maps` package via ggplot2)
## AoU `state` is a 2-letter abbreviation -> lowercase state name.
########################################################################

us_map <- map_data("state")

## Join per-disease state values onto the US polygons (one facet copy per
## disease). Suppressed / NA states keep their polygon row with val = NA.
poly_data <- function(tab, diseases, value_col) {
  d <- tab %>%
    filter(Disease %in% diseases) %>%
    mutate(map_region = tolower(state.name[match(state, state.abb)]),
           val = .data[[value_col]])
  dropped <- setdiff(unique(d$state[is.na(d$map_region)]), NA)
  if (length(dropped)) message("States not on the contiguous map (dropped): ",
                               paste(sort(dropped), collapse = ", "))
  us_map %>%
    rename(map_region = region) %>%
    inner_join(d %>% filter(!is.na(map_region)) %>% select(Disease, map_region, val),
               by = "map_region", relationship = "many-to-many") %>%
    mutate(Disease = factor(Disease, levels = diseases))
}

## Faceted choropleth: diverging fill (blue = below the midpoint, red = above),
## auto-ranged, centered at `mid`. NA / suppressed states are drawn with a
## diagonal HATCH (via ggpattern) instead of a flat grey.
choropleth_facet <- function(tab, diseases, value_col, mid, fill_lab, title) {
  dd   <- poly_data(tab, diseases, value_col)
  supp <- dplyr::filter(dd, is.na(val))

  p <- ggplot(dd, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = val), color = "grey55", linewidth = 0.08)

  if (nrow(supp) > 0) {
    if (requireNamespace("ggpattern", quietly = TRUE)) {
      p <- p + ggpattern::geom_polygon_pattern(
        data = supp, aes(long, lat, group = group),
        pattern = "stripe", pattern_angle = 45, pattern_density = 0.1,
        pattern_spacing = 0.015, pattern_size = 0.1,
        pattern_fill = "grey35", pattern_colour = NA,
        fill = "grey92", colour = "grey55", linewidth = 0.08)
    } else {                       # fallback if ggpattern is unavailable
      p <- p + geom_polygon(data = supp, aes(long, lat, group = group),
                            fill = "grey75", colour = "grey55", linewidth = 0.08)
    }
  }

  p +
    coord_quickmap() +
    facet_wrap(~ Disease, labeller = label_wrap_gen(width = 28)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = mid, na.value = "grey75", name = fill_lab) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(strip.text = element_text(size = 8, face = "bold"),
          legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5))
}

########################################################################
## Figures
########################################################################

dir.create("../graph", showWarnings = FALSE, recursive = TRUE)

## Each set gets three ratio maps: G-RAILS/unw, S-RAILS/unw, S-RAILS/G-RAILS
map_set <- function(diseases, tag, tag_label) {
  if (length(diseases) == 0) return(invisible(NULL))
  specs <- list(
    list(col = "r_g_unw", lab = "G-RAILS /\nunweighted", what = "G-RAILS vs Unweighted",
         file = paste0("map_", tag, "_global_vs_unw", OUTPUT_TAG, ".png")),
    list(col = "r_s_unw", lab = "S-RAILS /\nunweighted", what = "S-RAILS vs Unweighted",
         file = paste0("map_", tag, "_subgroup_vs_unw", OUTPUT_TAG, ".png")),
    list(col = "r_s_g",   lab = "S-RAILS /\nG-RAILS", what = "S-RAILS vs G-RAILS",
         file = paste0("map_", tag, "_subgroup_vs_global", OUTPUT_TAG, ".png"))
  )
  for (s in specs) {
    p <- choropleth_facet(state_tab, diseases, s$col, mid = 1, fill_lab = s$lab,
                          title = paste0(s$what, "  (", tag_label, ")"))
    print(p)
    ggsave(file.path("../graph", s$file), p, width = 12, height = 8, dpi = 300)
  }
}

map_set(setA,     "common", "most common")
map_set(setHighD, "highD",  "highest D")
map_set(setLowD,  "lowD",   "lowest D")
map_set(setHighE, "highE",  "highest E")
map_set(setLowE,  "lowE",   "lowest E")

########################################################################
## Regime example maps: one row per D-E regime, ordered high D / high E ->
## high D / low E -> low D / high E -> low D / low E, with that regime's
## N_EXAMPLES representative phenotypes (from the regime block above) across
## the row. Drawn for the three state-level ratios (S-RAILS/G-RAILS, the ratio
## D and E are built from; S-RAILS/unweighted and G-RAILS/unweighted as the
## parallel comparisons). The figure counterpart of the manuscript's
## representative-examples table. Each panel is titled with the phenotype (row 1)
## and the leverage subgroup (row 2): the state with the largest share of D, its
## share and its S/G ratio, prefixed "spread;" when the share is below
## LEVERAGE_SHARE_MIN.
## Colour: symmetric log scale centred at 1, squished to [1/MAP_RATIO_CAP,
## MAP_RATIO_CAP], so a single extreme state cannot flatten every other panel.
########################################################################

REGIME_ROW_ORDER <- c("High D, high E", "High D, low E", "Low D, high E", "Low D, low E")

regime_map <- function(value_col, fill_lab, what, file_tag) {
  rm_key <- regime_examples %>%
    mutate(regime = factor(as.character(regime), levels = REGIME_ROW_ORDER)) %>%
    group_by(regime) %>%
    mutate(slot = factor(paste0("Example ", row_number()),
                         levels = paste0("Example ", seq_len(N_EXAMPLES)))) %>%
    ungroup() %>%
    transmute(Disease, regime, slot,
              ## second title row: the leverage subgroup - the state with the
              ## largest share of D, its share and its S/G ratio, for every panel
              ## (e.g. "[CT ~95% D, S/G = 3.16]"); "spread;" is prefixed when the
              ## share is below LEVERAGE_SHARE_MIN (e.g. "[spread; NY ~33% D, S/G = 0.94]")
              state_label = ifelse(is.na(top_state), "",
                sprintf("[%s%s ~%.0f%% D, S/G = %.2f]",
                        ifelse(top_share >= LEVERAGE_SHARE_MIN, "", "spread; "),
                        top_state, 100 * top_share, top_state_ratio_SG)))

  rm_dd <- poly_data(state_tab, rm_key$Disease, value_col) %>%
    mutate(Disease = as.character(Disease)) %>%
    inner_join(rm_key, by = "Disease")
  rm_supp <- dplyr::filter(rm_dd, is.na(val))

  p <- ggplot(rm_dd, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = val), color = "grey55", linewidth = 0.08)
  if (nrow(rm_supp) > 0) {
    if (requireNamespace("ggpattern", quietly = TRUE)) {
      p <- p + ggpattern::geom_polygon_pattern(
        data = rm_supp, aes(long, lat, group = group),
        pattern = "stripe", pattern_angle = 45, pattern_density = 0.1,
        pattern_spacing = 0.015, pattern_size = 0.1,
        pattern_fill = "grey35", pattern_colour = NA,
        fill = "grey92", colour = "grey55", linewidth = 0.08)
    } else {
      p <- p + geom_polygon(data = rm_supp, aes(long, lat, group = group),
                            fill = "grey75", colour = "grey55", linewidth = 0.08)
    }
  }
  cap_breaks <- c(1 / MAP_RATIO_CAP, 0.67, 1, 1.5, MAP_RATIO_CAP)
  p +
    geom_text(data = rm_key, aes(x = -96, y = 55.5, label = Disease),
              inherit.aes = FALSE, size = 3.1, fontface = "bold", vjust = 1) +
    geom_text(data = rm_key, aes(x = -96, y = 53.2, label = state_label),
              inherit.aes = FALSE, size = 2.7, colour = "grey30", vjust = 1) +
    coord_quickmap(ylim = c(24, 56)) +
    facet_grid(regime ~ slot, switch = "y", drop = FALSE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, trans = "log2",
                         limits = c(1 / MAP_RATIO_CAP, MAP_RATIO_CAP),
                         oob = scales::squish, breaks = cap_breaks,
                         labels = c(paste0("≤", format(round(1 / MAP_RATIO_CAP, 2))),
                                    "0.67", "1", "1.5",
                                    paste0("≥", format(MAP_RATIO_CAP))),
                         na.value = "grey75", name = fill_lab) +
    labs(title = paste0(what, " state ratio by D-E regime"), x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(strip.text.x = element_text(size = 9, colour = "grey40"),
          strip.text.y.left = element_text(size = 10, face = "bold", angle = 90),
          strip.placement = "outside",
          legend.position = "right",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", hjust = 0.5))
}

if (nrow(regime_examples) > 0) {
  regime_specs <- list(
    list(col = "r_s_g",   lab = "S-RAILS /\nG-RAILS",    what = "S-RAILS vs G-RAILS",    tag = "SG"),
    list(col = "r_s_unw", lab = "S-RAILS /\nunweighted", what = "S-RAILS vs unweighted", tag = "SU"),
    list(col = "r_g_unw", lab = "G-RAILS /\nunweighted", what = "G-RAILS vs unweighted", tag = "GU")
  )
  for (s in regime_specs) {
    p_rm <- regime_map(s$col, s$lab, s$what, s$tag)
    print(p_rm)
    ggsave(paste0("../graph/map_DE_regimes_", s$tag, OUTPUT_TAG, ".png"), p_rm,
           width = 11, height = 10.5, dpi = 300)
  }
}

########################################################################
## Summary plots of the phenotype-level metrics (D_j, E_j)
########################################################################

## Highlight three sets with distinct colors (most common; highest D; highest E)
## and label only those; everything else is neutral grey. A phenotype in more
## than one extreme set is assigned once, in the order common > D > E.
LABEL_N <- 6
hiD <- pheno %>% filter(is.finite(D_js)) %>% slice_max(D_js, n = LABEL_N, with_ties = FALSE) %>% pull(Disease)
hiE <- pheno %>% filter(is.finite(E_js)) %>% slice_max(E_js, n = LABEL_N, with_ties = FALSE) %>% pull(Disease)
label_set <- union(hiD, hiE)

pheno_plot <- pheno %>%
  mutate(group = dplyr::case_when(
           Disease %in% setA ~ "most common",
           Disease %in% hiD  ~ "highest D",
           Disease %in% hiE  ~ "highest E",
           TRUE              ~ "other"),
         group = factor(group, levels = c("most common", "highest D", "highest E", "other")))

grp_cols <- c("most common" = "#1B7837", "highest D" = "#B2182B",
              "highest E" = "#2166AC", "other" = "grey80")

## (1) D vs E scatter — magnitude vs dispersion of the G-RAILS/S-RAILS difference
p_DE <- ggplot(pheno_plot, aes(D_js, E_js)) +
  geom_point(aes(color = group, size = P_unw), alpha = 0.85) +
  scale_color_manual(values = grp_cols,
                     breaks = c("most common", "highest D", "highest E"), name = NULL) +
  scale_size_continuous(name = "Unweighted\nprevalence", range = c(1.2, 7)) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4, alpha = 1)),
         size  = guide_legend(order = 2)) +
  labs(x = expression(D[j]), y = expression(E[j]),
       title = "Divergence D vs. dispersion E") +
  theme_bw(base_size = 13) +
  theme(legend.key = element_blank(), plot.title = element_text(face = "bold"))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_DE <- p_DE + ggrepel::geom_text_repel(
    data = dplyr::filter(pheno_plot, Disease %in% label_set),
    aes(label = Disease, color = group), size = 3, fontface = "plain",
    box.padding = 0.6, point.padding = 0.3, min.segment.length = 0,
    segment.size = 0.3, segment.color = "grey60", force = 2,
    max.overlaps = Inf, seed = 1, show.legend = FALSE)
}
print(p_DE)
ggsave(paste0("../graph/plot_D_vs_E", OUTPUT_TAG, ".png"), p_DE,
       width = 9.5, height = 7, dpi = 300)

## (2) Boxplots of D, log(D), E, and sum|log(S/G)| across phenotypes
box_df <- pheno %>%
  transmute(Disease,
            `D`             = D_js,
            `log(D)`        = log(D_js),
            `E`             = E_js,
            `sum|log(S/G)|` = sum_abs_logratio) %>%
  pivot_longer(-Disease, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("D", "log(D)", "E", "sum|log(S/G)|"))) %>%
  filter(is.finite(value))

p_box <- ggplot(box_df, aes(x = metric, y = value, fill = metric)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.4, alpha = 0.85) +
  geom_jitter(width = 0.12, alpha = 0.2, size = 0.6) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  facet_wrap(~ metric, scales = "free", nrow = 1) +
  labs(x = NULL, y = NULL, title = "Phenotype-level metric distributions") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(p_box)
ggsave(paste0("../graph/plot_D_E_boxplots", OUTPUT_TAG, ".png"), p_box,
       width = 11, height = 4.5, dpi = 300)

########################################################################
## Save tables (suppressed cells already NA in the ratio columns)
########################################################################

write_excel_csv(state_tab, paste0("../data/state_prevalence_ratios",   OUTPUT_TAG, ".csv"))
write_excel_csv(pheno,     paste0("../data/state_phenotype_selection", OUTPUT_TAG, ".csv"))
