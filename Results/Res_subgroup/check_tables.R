# check_tables.R -----------------------------------------------------------
# Reconcile the two simulation tables in the manuscript against the cached
# `fun.rep` summaries they were built from.
#
# Parses Table 1 (population, `tab1`) and Table 2 (subgroup, `tab2`) straight
# out of full_tex.tex, recomputes every cell from the history data, and writes
# table_vs_data_check.csv listing each cell with both values and the gap.
#
# Reported numbers are rounded to 2 dp in the paper, so anything under `tol`
# counts as agreement.

library(dplyr)

texfile <- "C:/Work/REpo/RAILS_Edit/Application/full_tex.tex"
datadir <- "C:/Work/2025 Fall/Subgroup/V12_new"
outdir  <- getwd()
tol     <- 0.011

cases   <- c("000", "100", "010", "001", "101", "111")   # Settings 1..6
labels  <- c("Naive", "Oracle", "G-Raking", "NPS", "GVS",
             "G-RAILS", "S-Raking", "S-RAILS")
groups  <- c("V1", "V2", "V3", "V4")

## ---- 1. parse the LaTeX tables ---------------------------------------------

tex <- readLines(texfile, warn = FALSE)

region <- function(start_pat, end_pat) {
  i <- grep(start_pat, tex, fixed = TRUE)[1]
  j <- grep(end_pat, tex[i:length(tex)], fixed = TRUE)[1] + i - 1
  paste(tex[i:j], collapse = " ")
}

# The setting number is not reliably recoverable from the row text (\cline and
# the "G1 Size" cells carry stray digits), so rows are indexed by position: the
# eight methods repeat in a fixed order, six times.
parse_rows <- function(txt, ncol) {
  rows <- strsplit(txt, "\\\\", fixed = TRUE)[[1]]      # split on the row break
  out <- NULL
  for (r in rows) {
    if (!grepl("textit", r)) next
    f <- strsplit(r, "&", fixed = TRUE)[[1]]
    k <- grep("textit", f)[1]
    meth <- sub(".*textit\\{([^}]*)\\}.*", "\\1", f[k])
    v <- f[-seq_len(k)]
    v <- gsub("cellcolor\\{[a-z]*\\}|cline\\{[0-9-]*\\}|[\\\\$~ ]", "", v)
    v <- sub("^<0\\.001$", "0.001", v)                  # "$< 0.001$" entries
    v <- suppressWarnings(as.numeric(v))
    v <- v[!is.na(v)]
    if (length(v) < ncol) next
    out <- rbind(out, data.frame(method = meth, t(v[seq_len(ncol)]),
                                 row.names = NULL))
  }
  out$setting <- rep(1:6, each = 8)[seq_len(nrow(out))]
  # the fixed method order is what makes the positional indexing safe
  stopifnot(identical(tolower(out$method), tolower(rep(labels, 6))))
  out[, c("setting", setdiff(names(out), c("setting", "method")), "method")][
    , c("setting", "method", setdiff(names(out), c("setting", "method")))]
}

t1 <- parse_rows(region("\\label{tab1}", "\\end{longtable}"), 6)
colnames(t1)[-(1:2)] <- c("RelBias", "AVar", "EVar", "MAD", "NomCP", "OraCP")

t2 <- parse_rows(region("\\label{tab2}", "\\end{tabular}"), 12)
colnames(t2)[-(1:2)] <- paste0(rep(c("RelBias", "EVar", "OrcCP"), 4),
                               "_G", rep(1:4, each = 3))

stopifnot(nrow(t1) == 48, nrow(t2) == 48)

## ---- 2. recompute the same quantities from the cached summaries -------------

d1 <- NULL; d2 <- NULL
for (i in seq_along(cases)) {
  o  <- readRDS(file.path(datadir, paste0("X5_V12_summary_", cases[i], ".RData")))
  pp <- as.numeric(read.csv(file.path(datadir, paste0("pre_pop_X5_V12_", cases[i], ".csv")))[, 2])
  m  <- o$out_mean

  d1 <- rbind(d1, data.frame(
    setting = i, method = labels,
    RelBias = m[["Relative Bias(Mean)"]] * 100, AVar = m[["Avg SD"]] * 1e3,
    EVar = m[["Emp SD"]] * 1e3, MAD = m[["MAD"]] * 1e3,
    NomCP = m[["TrueCoverage"]] * 100, OraCP = m[["Oracle Coverage"]] * 100,
    row.names = NULL))

  # Table 2's Rel Bias column reproduces what subgroup_summary.Rmd prints, NOT
  # bias/truth. The Rmd transposes first (line 476) and then sweeps over
  # MARGIN = 2, which is the *method* axis, so the length-4 pre_pop is recycled
  # across the eight method columns: each method is divided by the truth of
  # whichever subgroup its column position lands on. `rmd` replays that;
  # `correct` is bias divided by its own subgroup's truth.
  osb <- t(o$out_sub_bias); osb <- rbind(osb, Total = colSums(osb))
  rmd <- t(sweep(osb, 2, pp[-1], FUN = "/"))[, 1:4] * 1e4

  r <- data.frame(setting = i, method = labels, row.names = NULL)
  for (g in 1:4) {
    r[[paste0("RelBias_G", g)]] <- rmd[, g]                                      # as printed
    r[[paste0("EVar_G",    g)]] <- o$out_list$variance[, groups[g]] * 1e6
    r[[paste0("OrcCP_G",   g)]] <- o$out_sub_cov[, groups[g]] * 100
  }
  d2 <- rbind(d2, r)
}

## ---- 3. compare -------------------------------------------------------------

compare <- function(tab, dat, which_table) {
  cols <- setdiff(names(tab), c("setting", "method"))
  bind_rows(lapply(cols, function(cl)
    data.frame(table = which_table, setting = tab$setting, method = tab$method,
               column = cl, table_value = tab[[cl]],
               data_value = dat[[cl]][match(paste(tab$setting, tab$method),
                                            paste(dat$setting, dat$method))],
               row.names = NULL))) |>
    mutate(diff = table_value - data_value,
           agrees = abs(diff) < tol)
}

chk <- bind_rows(compare(t1, d1, "Table 1 (population)"),
                 compare(t2, d2, "Table 2 (subgroup)"))
write.csv(chk, file.path(outdir, "table_vs_data_check.csv"), row.names = FALSE)

cat("\n=== agreement by table and column ===\n")
print(chk |> group_by(table, column) |>
        summarise(n = n(), agree = sum(agrees), mismatch = sum(!agrees),
                  max_abs_diff = round(max(abs(diff)), 2), .groups = "drop") |>
        as.data.frame())

cat("\n=== the 20 largest mismatches ===\n")
print(chk |> filter(!agrees) |> arrange(desc(abs(diff))) |> head(20) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))) |> as.data.frame())

cat("\nTable 2 RelBias columns above are on the 1e-4 scale used by the table.",
    "\nOn a true percentage scale the same cells are 100x smaller",
    "(e.g. Setting 1 / Naive / G1 = 4.59%, printed as 459.62).\n")

cat("\n=== Table 1: every mismatching cell ===\n")
print(chk |> filter(table == "Table 1 (population)", !agrees) |>
        arrange(setting, method) |>
        mutate(across(where(is.numeric), \(x) round(x, 3))) |> as.data.frame())

cat("\n=== Table 2: agreement by setting ===\n")
print(chk |> filter(table == "Table 2 (subgroup)") |>
        mutate(kind = sub("_G[0-9]$", "", column)) |>
        group_by(setting, kind) |>
        summarise(agree = sum(agrees), n = n(), .groups = "drop") |>
        tidyr::pivot_wider(names_from = kind, values_from = c(agree, n)) |>
        as.data.frame())

# Are whole blocks of Table 2 duplicated between settings? Compare each pair of
# settings cell by cell, using the printed values only.
cat("\n=== Table 2: identical cells between pairs of settings (out of 96) ===\n")
num <- setdiff(names(t2), c("setting", "method"))
dup <- expand.grid(a = 1:6, b = 1:6) |> filter(a < b)
dup$identical_cells <- mapply(function(a, b)
  sum(as.matrix(t2[t2$setting == a, num]) == as.matrix(t2[t2$setting == b, num]),
      na.rm = TRUE), dup$a, dup$b)
print(dup[order(-dup$identical_cells), ][1:6, ], row.names = FALSE)

cat("\nwrote ", file.path(outdir, "table_vs_data_check.csv"), "\n", sep = "")
