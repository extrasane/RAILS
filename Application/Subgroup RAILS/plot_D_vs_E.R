#### plot_D_vs_E.R — the D_j vs E_j phenotype scatter
####
#### Drop-in replacement for the summary-scatter block of
#### 10_State_Prevalence_Maps.R. Source this file, then call
####
####   p_DE <- plot_D_vs_E(pheno, setA, setHighD, setLowD, setHighE, setLowE)
####   ggsave("../graph/plot_D_vs_E.png", p_DE, width = 9, height = 6.2, dpi = 300)
####
#### x = D_j, the population-weighted G-vs-S Jensen-Shannon divergence.
#### y = E_j, its dispersion entropy across states (1 = spread evenly over all
####     states, 0 = concentrated in a few).
####
#### Colour marks the four extreme sets that get their own map figures, so the
#### reader can find each mapped phenotype in the scatter. A phenotype can fall
#### in more than one set; membership is resolved by the priority below and the
#### number of multiply-assigned phenotypes is reported, rather than silently
#### picking one. "Most common" is drawn as a ring instead of a fill, so it can
#### coexist with a D/E extreme on the same point.

library(ggplot2)
library(dplyr)

SET_COLOURS <- c(
  "highest D" = "#C0392B",   # largest G-vs-S divergence
  "lowest D"  = "#2C7FB8",   # smallest
  "highest E" = "#1B7837",   # divergence spread evenly across states
  "lowest E"  = "#E08214",   # divergence concentrated in a few states
  "other"     = "grey82"
)

plot_D_vs_E <- function(pheno, setA, setHighD, setLowD, setHighE, setLowE,
                        n_label_D = 5, n_label_E_hi = 3, n_label_E_lo = 3,
                        wrap_width = 24, verbose = TRUE) {

  ## ---- set membership, priority-resolved ---------------------------------
  d <- pheno %>%
    mutate(
      set = case_when(
        Disease %in% setHighD ~ "highest D",
        Disease %in% setLowD  ~ "lowest D",
        Disease %in% setHighE ~ "highest E",
        Disease %in% setLowE  ~ "lowest E",
        TRUE                  ~ "other"),
      set    = factor(set, levels = names(SET_COLOURS)),
      common = ifelse(Disease %in% setA, "most common", "other")
    )

  if (verbose) {
    n_multi <- sum(rowSums(cbind(d$Disease %in% setHighD, d$Disease %in% setLowD,
                                 d$Disease %in% setHighE, d$Disease %in% setLowE)) > 1)
    message("D/E extreme sets: ", sum(d$set != "other"), " phenotypes coloured, ",
            n_multi, " of them in more than one set (shown by the priority ",
            "highest D > lowest D > highest E > lowest E).")
  }

  ## ---- what to label ------------------------------------------------------
  ## Labelling the top E_j outright picks a pile of points stacked on the D = 0
  ## edge, which is what made the old version collide. Take the high-E labels
  ## from the phenotypes that also have some divergence to speak of.
  med_D <- median(d$D_js, na.rm = TRUE)
  lab <- bind_rows(
    d %>% filter(is.finite(D_js)) %>% slice_max(D_js, n = n_label_D, with_ties = FALSE),
    d %>% filter(is.finite(E_js), D_js >= med_D) %>%
      slice_max(E_js, n = n_label_E_hi, with_ties = FALSE),
    d %>% filter(is.finite(E_js)) %>% slice_min(E_js, n = n_label_E_lo, with_ties = FALSE)
  ) %>% distinct(Disease, .keep_all = TRUE)

  ## ---- plot ---------------------------------------------------------------
  p <- ggplot(d, aes(D_js, E_js)) +
    geom_point(aes(fill = set, size = P_unw, colour = common),
               shape = 21, alpha = 0.85, stroke = 0.6) +
    scale_fill_manual(values = SET_COLOURS, name = "phenotype set",
                      breaks = setdiff(names(SET_COLOURS), "other")) +
    scale_colour_manual(values = c("most common" = "grey15", "other" = NA),
                        breaks = "most common", name = NULL,
                        na.value = NA) +
    scale_size_continuous(name = "unweighted\nprevalence", range = c(1.6, 7),
                          breaks = c(0.05, 0.15, 0.30),
                          labels = scales::percent_format(accuracy = 1)) +
    ## 1e-5 units on the axis instead of 0e+00 / 2e-05 tick labels
    scale_x_continuous(labels = function(v) format(v * 1e5, trim = TRUE),
                       expand = expansion(mult = c(0.06, 0.10))) +
    labs(x = expression(D[j] ~ " (" %*% 10^-5 * ")"), y = expression(E[j]),
         title = "Divergence and dispersion of the G-RAILS / S-RAILS difference",
         subtitle = paste("Right = the two weightings disagree more;",
                          "low = that disagreement sits in a few states")) +
    guides(
      fill   = guide_legend(order = 1, override.aes = list(size = 4, colour = NA)),
      colour = guide_legend(order = 2, override.aes = list(size = 4, fill = "grey82",
                                                           stroke = 0.9)),
      size   = guide_legend(order = 3, override.aes = list(fill = "grey70", colour = NA))
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, colour = "grey35"),
          legend.key = element_blank(),
          legend.spacing.y = unit(2, "pt"))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = lab, aes(label = stringr::str_wrap(Disease, wrap_width)),
      size = 3.1, lineheight = 0.92, colour = "grey15",
      segment.colour = "grey55", segment.size = 0.3, min.segment.length = 0,
      box.padding = 0.55, point.padding = 0.35, force = 3,
      max.overlaps = Inf, seed = 1, show.legend = FALSE)
  }
  p
}
