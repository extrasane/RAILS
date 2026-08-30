library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
scenario <- "X5_V12"
cases <- c("000","001","010","011","100","101","110","111")
level_order <- c("V1", "V2", "V3", "V4")
load_one_case <- function(case,
                          scenario,
                          sim_result,
                          n,
                          pre_pop,
                          truey,
                          var_method,
                          var_stratified) {
  
  fname <- paste0(scenario, "_summary_", case, ".RData")
  
  if (fname %in% list.files(getwd(), recursive = TRUE)) {
    out <- readRDS(fname)
  } else {
    out <- fun.rep(
      x = sim_result,
      n = n,
      pre_pop = pre_pop,
      truey = truey,
      var_method = var_method,
      var_stratified = var_stratified
    )
    saveRDS(out, fname)
  }
  
  out$out_list$Coverage
}

coverage_list <- map(
  cases,
  load_one_case,
  scenario = scenario,
  sim_result = sim_result,
  n = n,
  pre_pop = pre_pop,
  truey = truey,
  var_method = var_method,
  var_stratified = var_stratified
)

names(coverage_list) <- cases
coverage_list <- map(
  coverage_list,
  ~ .x[, level_order, drop = FALSE]
)



coverage_long_all <- imap_dfr(
  coverage_list,
  function(df, case) {
    
    df %>%
      as.data.frame() %>%
      rownames_to_column("Method") %>%
      pivot_longer(
        cols = all_of(level_order),   # ← critical fix
        names_to = "V",
        values_to = "Coverage"
      ) %>%
      mutate(
        Case = case,
        Case_V = paste(case, V, sep = "_"),
      )
  }
)
case_v_levels <- as.vector(
  outer(cases, level_order, paste, sep = "_")
)
case_v_method <- as.vector(
  outer(cases, rownames(coverage_list[[1]]), paste, sep = "_")
)


coverage_long_all <- coverage_long_all %>%
  mutate(
    Case_V = factor(Case_V, levels = case_v_levels),
    Method = factor(Method, levels = rownames(coverage_list[[1]]))
  )

ggplot(coverage_long_all,
       aes(x = Case_V, y = Method, fill = Coverage)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_y_discrete(limits = rev(rownames(coverage_list[[1]]))) +  # Add this line
  scale_fill_gradient2(
    low = "firebrick",     
    mid = "white",
    high = "darkblue",
    midpoint = 0.95,
    limits = c(0, 1),
    name = "Coverage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Scenario by Variable",
    y = "Method"
  )

