# Validate the oracle-coverage routine by applying the SAME code to the
# population, where fun.rep already cached an answer to compare against.

datadir <- "C:/Work/2025 Fall/Subgroup/V12_new"
cases   <- c("000", "100", "010", "001", "101", "111")
labels  <- c("Naive","Oracle","G-Raking","NPS","GVS","G-RAILS","S-Raking","S-RAILS")
raw_key <- c("unweighted","trueweight","svy_oneway_rake","nps","GVS",
             "rails","rake_subgroup","rails_subgroup")

for (i in seq_along(cases)) {
  o  <- readRDS(file.path(datadir, paste0("X5_V12_summary_", cases[i], ".RData")))
  pp <- as.numeric(read.csv(file.path(datadir,
          paste0("pre_pop_X5_V12_", cases[i], ".csv")))[, 2])
  truey <- pp[1]
  sr <- readRDS(file.path(datadir, paste0("X5_V12_", cases[i], ".RData")))

  mb  <- o$out_mean[["Bias(Mean - True)"]]          # avg bias per method, population
  num <- den <- numeric(length(labels))
  for (rp in sr) {
    for (j in seq_along(labels)) {
      p <- rp[[raw_key[j]]]$pop
      if (is.null(p) || is.na(p[1]) || is.na(p[2])) next
      est <- p[1] - mb[j]                            # remove that method's avg bias
      ci  <- est + c(-1, 1) * 1.96 * sqrt(p[2])      # rebuild the interval
      num[j] <- num[j] + as.numeric(truey >= ci[1] & truey <= ci[2])
      den[j] <- den[j] + 1
    }
  }
  mine   <- num / den * 100
  cached <- o$out_mean[["Oracle Coverage"]] * 100
  cat(sprintf("Setting %d (%s)  max |mine - cached| = %.4f pp   exact matches %d/8\n",
              i, cases[i], max(abs(mine - cached)), sum(abs(mine - cached) < 1e-9)))
  if (i == 1) {
    cmp <- data.frame(method = labels, mine = round(mine, 2), cached = round(cached, 2))
    print(cmp, row.names = FALSE)
  }
}
