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
         file = paste0("map_", tag, "_global_vs_unw.png")),
    list(col = "r_s_unw", lab = "S-RAILS /\nunweighted", what = "S-RAILS vs Unweighted",
         file = paste0("map_", tag, "_subgroup_vs_unw.png")),
    list(col = "r_s_g",   lab = "S-RAILS /\nG-RAILS", what = "S-RAILS vs G-RAILS",
         file = paste0("map_", tag, "_subgroup_vs_global.png"))
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
ggsave("../graph/plot_D_vs_E.png", p_DE, width = 9.5, height = 7, dpi = 300)

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
ggsave("../graph/plot_D_E_boxplots.png", p_box, width = 11, height = 4.5, dpi = 300)

########################################################################
## Save tables (suppressed cells already NA in the ratio columns)
########################################################################

write_excel_csv(state_tab, "../data/state_prevalence_ratios.csv")
write_excel_csv(pheno,     "../data/state_phenotype_selection.csv")
