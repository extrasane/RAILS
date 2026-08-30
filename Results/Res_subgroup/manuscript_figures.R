# manuscript_figures.R -----------------------------------------------------
# Figures accompanying the two simulation tables in the subgroup manuscript.
#
#   fig_sim_population.png       companion to Table 1 (population estimation)
#   fig_sim_subgroup.png         companion to Table 2 (subgroup estimation)
#   fig_sim_bias_population.png  per-replicate bias, violin + box, population
#   fig_sim_bias_subgroup.png    per-replicate bias, violin + box, by subgroup
#
# Reads the cached `fun.rep` summaries (X5_V12_summary_<case>.RData) and the
# truth CSVs written by the run side. The subgroup oracle coverage is the one
# quantity `fun.rep` does not cache, so it is recomputed once from the raw
# replicate files (X5_V12_<case>.RData) and cached in oracle_cov_subgroup.csv.
#
# Setting numbering follows the manuscript: the six reported settings are the
# 3-bit cases below, in this order.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

datadir <- "C:/Work/2025 Fall/Subgroup/V12_new"   # where the cached results live
outdir  <- getwd()

cases   <- c("000", "100", "010", "001", "101", "111")   # -> Settings 1..6
groups  <- c("V1", "V2", "V3", "V4")

# `labels` is the order the cached matrices are stored in - do not reorder it,
# it is used positionally. `display` is the order the figures show.
methods <- c("naive", "oracle", "G-Raking", "NPS", "GVS",
             "G-RAILS", "S-Raking", "S-RAILS")
labels  <- c("Naive", "Oracle", "G-Raking", "NPS", "GVS",
             "G-RAILS", "S-Raking", "S-RAILS")
raw_key <- c("unweighted", "trueweight", "svy_oneway_rake", "nps", "GVS",
             "rails", "rake_subgroup", "rails_subgroup")
display <- c("Naive", "NPS", "GVS", "G-Raking", "S-Raking",
             "G-RAILS", "S-RAILS", "Oracle")

## ---- 1. subgroup oracle coverage (recomputed once, then cached) -------------
# Same construction fun.rep uses for the population: recentre each replicate's
# estimate by the method's average bias, rebuild the interval from that
# replicate's variance, and ask whether it covers the true subgroup prevalence.
# It isolates the variance estimate by taking the bias out.

oracle_file <- file.path(outdir, "oracle_cov_subgroup.csv")

if (file.exists(oracle_file)) {
  oracle_cov <- read.csv(oracle_file)
} else {
  message("computing subgroup oracle coverage from the raw replicate files ...")
  oracle_cov <- NULL
  for (i in seq_along(cases)) {
    o  <- readRDS(file.path(datadir, paste0("X5_V12_summary_", cases[i], ".RData")))
    pp <- as.numeric(read.csv(file.path(datadir,
            paste0("pre_pop_X5_V12_", cases[i], ".csv")))[, 2])[-1]
    mb <- as.matrix(o$out_sub_bias[, groups])          # mean bias, methods x groups
    sr <- readRDS(file.path(datadir, paste0("X5_V12_", cases[i], ".RData")))

    num <- den <- matrix(0, length(labels), length(groups))
    for (rp in sr) {
      for (j in seq_along(labels)) {
        s <- rp[[raw_key[j]]]$sub
        if (is.null(s) || !all(c("mean", "variance") %in% rownames(s))) next
        m <- as.numeric(s["mean", groups]); v <- as.numeric(s["variance", groups])
        ok <- !is.na(m) & !is.na(v)
        if (!any(ok)) next
        est <- m - mb[j, ]
        hit <- as.numeric(pp >= est - 1.96 * sqrt(v) & pp <= est + 1.96 * sqrt(v))
        num[j, ok] <- num[j, ok] + hit[ok]
        den[j, ok] <- den[j, ok] + 1
      }
    }
    oracle_cov <- rbind(oracle_cov, data.frame(
      setting = i, method = rep(labels, length(groups)),
      group = rep(paste0("G", seq_along(groups)), each = length(labels)),
      OraCP = as.vector(num / den * 100), row.names = NULL))
  }
  write.csv(oracle_cov, oracle_file, row.names = FALSE)
}

## ---- 2. gather --------------------------------------------------------------

pop <- NULL; sub <- NULL; rep_pop <- NULL; rep_sub <- NULL

for (i in seq_along(cases)) {
  cs <- cases[i]
  o  <- readRDS(file.path(datadir, paste0("X5_V12_summary_", cs, ".RData")))
  pp <- as.numeric(read.csv(file.path(datadir, paste0("pre_pop_X5_V12_", cs, ".csv")))[, 2])
  truey <- pp[1]; truesub <- pp[-1]

  rep_pop <- rbind(rep_pop, data.frame(setting = i, case = cs,
                                       method = labels[match(o$out_p_all$Method, methods)],
                                       value  = o$out_p_all$value, row.names = NULL))
  for (g in seq_along(groups)) {
    pg <- o$out_p_sub_bias[[groups[g]]]
    rep_sub <- rbind(rep_sub, data.frame(setting = i, case = cs, group = paste0("G", g),
                                         method = labels[match(pg$Method, methods)],
                                         value  = pg$value, row.names = NULL))
  }

  m <- o$out_mean
  pop <- rbind(pop, data.frame(
    setting  = i, case = cs, method = labels,
    rel_bias = m[["Relative Bias(Mean)"]] * 100,          # per cent
    AVar     = m[["Avg SD"]] * 1e3,
    EVar     = m[["Emp SD"]] * 1e3,
    MAD      = m[["MAD"]]    * 1e3,
    NomCP    = m[["TrueCoverage"]]    * 100,
    OraCP    = m[["Oracle Coverage"]] * 100,
    truth    = truey, row.names = NULL))

  for (g in seq_along(groups)) {
    # note: bias is divided by its OWN subgroup's truth (see check_tables.R for
    # why Table 2's column does not match this)
    sub <- rbind(sub, data.frame(
      setting  = i, case = cs, method = labels, group = paste0("G", g),
      rel_bias = o$out_sub_bias[, groups[g]] / truesub[g] * 100,
      EVar     = o$out_list$variance[, groups[g]] * 1e6,
      NomCP    = o$out_sub_cov[, groups[g]] * 100,
      size     = o$out_list$size[, groups[g]],
      truth    = truesub[g], row.names = NULL))
  }
}

sub <- left_join(sub, oracle_cov, by = c("setting", "method", "group"))

fct <- function(d) transform(d,
  method  = factor(method, levels = display),
  setting = factor(setting, levels = 1:6,
                   labels = paste0("Setting ", 1:6, "\n(", cases, ")")))
pop <- fct(pop); sub <- fct(sub); rep_pop <- fct(rep_pop); rep_sub <- fct(rep_sub)

write.csv(pop, file.path(outdir, "tidy_pop.csv"), row.names = FALSE)
write.csv(sub, file.path(outdir, "tidy_sub.csv"), row.names = FALSE)

## ---- 3. shared style --------------------------------------------------------

# Naive / Oracle are references (grey, black); global methods blue;
# subgroup methods orange-red.
pal <- c("Naive"    = "#9aa5ae", "Oracle"   = "#16212b",
         "G-Raking" = "#9ecae1", "NPS"      = "#6baed6",
         "GVS"      = "#3182bd", "G-RAILS"  = "#08519c",
         "S-Raking" = "#fdae6b", "S-RAILS"  = "#c4451c")
shp <- c("Naive" = 4, "Oracle" = 23, "G-Raking" = 21, "NPS" = 21,
         "GVS" = 21, "G-RAILS" = 21, "S-Raking" = 21, "S-RAILS" = 21)

base <- theme_bw(base_size = 11, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "#eef2f7", colour = "#c3ccd6"),
        strip.text = element_text(size = 9, lineheight = 1.05),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(t = -2),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.6, colour = "#5d6b79"))

dodge <- position_dodge(width = 0.78)
sep   <- geom_vline(xintercept = seq(1.5, 5.5, 1), colour = "#dde3ea", linewidth = .4)

## ---- 4. Figure 1: population ------------------------------------------------

p1a <- ggplot(pop, aes(setting, pmax(abs(rel_bias), 0.01),
                       fill = method, shape = method, colour = method)) +
  sep +
  geom_point(position = dodge, size = 2.1, stroke = .5) +
  scale_y_log10(breaks = c(.01, .03, .1, .3, 1, 3, 10),
                labels = c("0.01", "0.03", "0.1", "0.3", "1", "3", "10")) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_shape_manual(values = shp) +
  labs(y = "|relative bias|  (%)", x = NULL,
       title = "(a)  Population prevalence — absolute relative bias",
       subtitle = "log scale, lower is better; values under 0.01% are drawn at 0.01") +
  base

p1b <- ggplot(pop, aes(setting, NomCP, fill = method, shape = method, colour = method)) +
  sep +
  geom_hline(yintercept = 95, linetype = "22", colour = "#c4451c", linewidth = .45) +
  geom_point(position = dodge, size = 2.1, stroke = .5) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_shape_manual(values = shp) +
  labs(y = "nominal coverage  (%)", x = NULL,
       title = "(b)  Population prevalence — nominal coverage",
       subtitle = "dashed line: 95% nominal level") +
  base

fig1 <- p1a / p1b + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "fig_sim_population.png"), fig1,
       width = 10, height = 7, dpi = 300, bg = "white")

## ---- 5. Figure 2: subgroup --------------------------------------------------
# Rows run Naive .. Oracle top to bottom, with the two reference methods ruled
# off from the six competing ones.

heat_base <- base + theme(panel.grid = element_blank(),
                          legend.position = "right",
                          legend.title = element_text(size = 8.5))
rule <- geom_hline(yintercept = c(1.5, 7.5), colour = "#16212b", linewidth = .5)

# |relative bias|: Naive runs to ~12% and would flatten everything else, so the
# fill is capped at 4% with out-of-range cells squished to the top of the ramp.
p2a <- ggplot(sub, aes(group, method, fill = pmin(abs(rel_bias), 4))) +
  geom_tile(colour = "white", linewidth = .35) + rule +
  facet_grid(~ setting) +
  scale_y_discrete(limits = rev(display)) +
  scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 4),
                      oob = scales::squish, name = "|rel. bias|\n(%)") +
  labs(x = NULL, y = NULL,
       title = "(a)  Subgroup prevalence — absolute relative bias",
       subtitle = "white = unbiased, red = biased; scale capped at 4%") +
  heat_base

cov_panel <- function(v, ttl, sub_ttl, xlab = NULL) {
  ggplot(sub, aes(group, method, fill = .data[[v]])) +
    geom_tile(colour = "white", linewidth = .35) + rule +
    facet_grid(~ setting) +
    scale_y_discrete(limits = rev(display)) +
    scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 100),
                        name = "coverage\n(%)") +
    labs(x = xlab, y = NULL, title = ttl, subtitle = sub_ttl) +
    heat_base
}

p2b <- cov_panel("NomCP", "(b)  Subgroup prevalence — nominal coverage",
                 "covers the true subgroup prevalence; deep red = at the 95% level, white = no coverage")
p2c <- cov_panel("OraCP", "(c)  Subgroup prevalence — oracle coverage",
                 "the same interval recentred by the method's average bias, so it reflects the variance alone",
                 xlab = "subgroup")

fig2 <- p2a / p2b / p2c
ggsave(file.path(outdir, "fig_sim_subgroup.png"), fig2,
       width = 11, height = 10.5, dpi = 300, bg = "white")

## ---- 6. Figures 3-4: per-replicate bias, violin + box -----------------------
# Same idiom as the main-paper simulation figure: a violin of the 1000
# replicate biases with a narrow box inside, methods on the vertical axis.

box_theme <- theme_bw(base_size = 13, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef2f7", colour = "#c3ccd6"),
        strip.text = element_text(size = 10.5, lineheight = 1.05),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9.5, colour = "#5d6b79"))

bias_panel <- function(d, ylim) {
  list(plot = ggplot(d, aes(method, value)) +
     geom_violin(aes(fill = method), width = 1, colour = NA, alpha = .6,
                 trim = FALSE, position = position_dodge(.9)) +
     geom_boxplot(width = .15, alpha = .2, outlier.size = .35,
                  outlier.alpha = .35, position = position_dodge(.9)) +
     geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
     scale_fill_manual(values = pal) +
     scale_x_discrete(limits = rev(display)) +
     coord_flip(ylim = ylim) +
     labs(x = "Method", y = "Bias"),
     drop = sum(d$value < ylim[1] | d$value > ylim[2]))
}

b3 <- bias_panel(rep_pop, c(-0.020, 0.035))
fig3 <- b3$plot + facet_wrap(~ setting, nrow = 1) + box_theme +
  labs(title = "Population prevalence — distribution of the bias over 1000 replicates",
       subtitle = sprintf("dashed line: zero bias; %d of %d replicate values fall outside the axis and are not drawn",
                          b3$drop, nrow(rep_pop)))
ggsave(file.path(outdir, "fig_sim_bias_population.png"), fig3,
       width = 34, height = 13, units = "cm", dpi = 300, bg = "white")

b4 <- bias_panel(rep_sub, c(-0.035, 0.045))
fig4 <- b4$plot + facet_grid(group ~ setting) + box_theme +
  labs(title = "Subgroup prevalence — distribution of the bias over 1000 replicates",
       subtitle = sprintf("rows are the subgroups defined by the separator; %d of %d replicate values fall outside the axis and are not drawn",
                          b4$drop, nrow(rep_sub)))
ggsave(file.path(outdir, "fig_sim_bias_subgroup.png"), fig4,
       width = 34, height = 30, units = "cm", dpi = 300, bg = "white")

message("wrote figures and tidy CSVs to ", normalizePath(outdir))
