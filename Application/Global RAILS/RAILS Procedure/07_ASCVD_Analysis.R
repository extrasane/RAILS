#### 07_ASCVD_Analysis.R — ASCVD risk distribution, unweighted vs RAILS
#### Run inside the AoU Researcher Workbench.
####
#### Compares the distribution of 10-year ASCVD risk before and after RAILS
#### weighting, and saves publication-quality figures.
####
#### Inputs:
####   raking_weights_*.csv — person IDs + weight columns (needs w_rails)
####   ASCVD_risks.csv      — person IDs + `risk`

library(tidyverse)

########################################################################
## Load data
########################################################################

## Bucket route (uncomment in the Workbench):
# my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
# for (f in c("raking_weights_v7.csv", "ASCVD_risks.csv")) {
#   system(paste0("gsutil cp ", my_bucket, "/data/", f, " ."), intern = TRUE)
# }
# raking.wts  <- read_csv("raking_weights_v7.csv")
# ASCVD_risks <- read_csv("ASCVD_risks.csv")

raking.wts  <- read_csv("../Andrew/jasa_raking_wts_w_diseases_2_code_requirement.csv")
ASCVD_risks <- read_csv("../data/ASCVD_risks.csv")

data <- ASCVD_risks %>%
  inner_join(raking.wts, by = "person_id")

## Long form: one row per person per weighting scheme
data2 <- data %>%
  select(person_id, risk, w_rails) %>%
  mutate(wt.class = "RAILS-Weighted") %>%
  rename(wt = w_rails) %>%
  rbind(data %>%
          select(person_id, risk) %>%
          mutate(wt.class = "Unweighted", wt = 1))

########################################################################
## Risk-category comparison — the numbers quoted in the manuscript
##   Risk categories follow the usual ASCVD convention:
##     low <5%, borderline/intermediate 5-10%, high >10%
##   Weighted proportion = sum(w * I(category)) / sum(w)
########################################################################

wtd_prop <- function(ind, w) sum(w * ind) / sum(w)

risk_summary <- tibble::tibble(
  Category = c("Low (<5%)", "Intermediate (5-10%)", "High (>10%)"),
  Unweighted = c(
    wtd_prop(data$risk < 5,                        rep(1, nrow(data))),
    wtd_prop(data$risk >= 5 & data$risk <= 10,     rep(1, nrow(data))),
    wtd_prop(data$risk > 10,                       rep(1, nrow(data)))
  ),
  RAILS = c(
    wtd_prop(data$risk < 5,                    data$w_rails),
    wtd_prop(data$risk >= 5 & data$risk <= 10, data$w_rails),
    wtd_prop(data$risk > 10,                   data$w_rails)
  )
) %>%
  mutate(
    Unweighted_pct = 100 * Unweighted,
    RAILS_pct      = 100 * RAILS,
    ## Percentage-POINT change (not a relative percent change)
    Diff_pp        = RAILS_pct - Unweighted_pct,
    ## Relative change, for contrast
    Rel_pct_change = 100 * (RAILS - Unweighted) / Unweighted
  )

print(risk_summary %>%
        select(Category, Unweighted_pct, RAILS_pct, Diff_pp, Rel_pct_change) %>%
        mutate(across(where(is.numeric), ~round(.x, 2))))

## Mean and SD of risk under each weighting
risk_moments <- tibble::tibble(
  Statistic  = c("Mean risk (%)", "SD of risk"),
  Unweighted = c(mean(data$risk), sd(data$risk)),
  RAILS      = c(weighted.mean(data$risk, data$w_rails),
                 sqrt(Hmisc::wtd.var(data$risk, data$w_rails)))
) %>%
  mutate(Diff = RAILS - Unweighted)
print(risk_moments %>% mutate(across(where(is.numeric), ~round(.x, 3))))

## Sentence-ready values
high_pp <- risk_summary$Diff_pp[risk_summary$Category == "High (>10%)"]
low_pp  <- risk_summary$Diff_pp[risk_summary$Category == "Low (<5%)"]
message(sprintf(
  "High-risk (>10%%): %.1f%% -> %.1f%% (%+.1f percentage points)\nLow-risk (<5%%):  %.1f%% -> %.1f%% (%+.1f percentage points)",
  risk_summary$Unweighted_pct[3], risk_summary$RAILS_pct[3], high_pp,
  risk_summary$Unweighted_pct[1], risk_summary$RAILS_pct[1], low_pp))

########################################################################
## Figure export helper — publication-quality TIFF + vector PDF + PNG.
## Text size, not dpi, is what makes manuscript figures look blurry, so the
## base font is set explicitly before saving.
########################################################################

save_pub <- function(plot, stem, width = 7, height = 5, base_size = 11) {
  p <- plot + theme(text        = element_text(size = base_size),
                    plot.title  = element_text(hjust = 0.5, size = base_size + 1),
                    axis.title  = element_text(size = base_size),
                    axis.text   = element_text(size = base_size - 2),
                    legend.text = element_text(size = base_size - 1))

  ggsave(paste0(stem, ".tiff"), plot = p, width = width, height = height,
         units = "in", dpi = 600, device = "tiff", compression = "lzw")
  ggsave(paste0(stem, ".pdf"),  plot = p, width = width, height = height,
         units = "in", device = cairo_pdf)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(paste0(stem, ".png"), plot = p, width = width, height = height,
           units = "in", dpi = 600, device = ragg::agg_png)
  } else {
    ggsave(paste0(stem, ".png"), plot = p, width = width, height = height,
           units = "in", dpi = 600, type = "cairo")
  }
  invisible(p)
}

########################################################################
## (1) Unweighted risk distribution
########################################################################

p_unw <- ggplot(data, aes(x = risk)) +
  geom_histogram(aes(y = after_stat(density)), color = "black", fill = "royalblue") +
  scale_x_continuous(trans = "log10") +
  labs(x = "Unweighted Risk (log10)", y = "Density",
       title = "Distribution of Unweighted Risk (log10-transformed)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(p_unw)
save_pub(p_unw, "ASCVD_Fig1_unweighted")

########################################################################
## (2) RAILS-weighted risk distribution
########################################################################

p_w <- ggplot(data, aes(x = risk)) +
  geom_histogram(aes(y = after_stat(density), weight = w_rails),
                 color = "black", fill = "royalblue") +
  scale_x_continuous(trans = "log10") +
  labs(x = "Weighted Risk (log10)", y = "Density",
       title = "Distribution of RAILS-Weighted Risk (log10-transformed)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(p_w)
save_pub(p_w, "ASCVD_Fig2_rails_weighted")

########################################################################
## (3) Overlaid histograms + density — unweighted vs RAILS
########################################################################

p_overlay <- ggplot(data2, aes(x = risk, fill = wt.class)) +
  geom_histogram(aes(y = after_stat(density), weight = wt),
                 color = "black", position = "identity", alpha = 0.5) +
  geom_density(aes(weight = wt), alpha = 0.25) +
  scale_x_continuous(trans = "log10") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "ASCVD Risk (log10)", y = "Density", fill = "Class",
       title = "Distribution of ASCVD Risk (log10-transformed)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(p_overlay)
save_pub(p_overlay, "ASCVD_Fig3_overlay_hist")

########################################################################
## (4) Density only — the cleanest version for the manuscript
########################################################################

p_density <- ggplot(data2, aes(x = risk, fill = wt.class)) +
  geom_density(aes(weight = wt), alpha = 0.5) +
  scale_x_continuous(trans = "log10") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "ASCVD Risk (log10)", y = "Density", fill = "Class",
       title = "Distribution of ASCVD Risk (log10-transformed)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(p_density)
save_pub(p_density, "ASCVD_Fig4_density")

########################################################################
## Optionally copy figures to the bucket
########################################################################

# system(paste0("gsutil cp ./ASCVD_Fig*.* ", Sys.getenv("WORKSPACE_BUCKET"), "/data/"), intern = TRUE)
