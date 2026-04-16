spec_color2 <- function(x, alpha = 1, begin = 0, end = 1,
                        direction = 1, option = "D",
                        na_color = "#BBBBBB", scale_from = NULL,
                        palette = viridisLite::viridis(256, alpha, begin, end, direction, option)) {
  n <- length(palette)
  if (is.null(scale_from)) {
    x <- round(scales::rescale(x, c(1, n)))
  } else {
    x <- round(scales::rescale(x, to = c(1, n),
                               from = scale_from))
  }
  
  color_code <- palette[x]
  color_code[is.na(color_code)] <- na_color
  return(color_code)
}

# Function to apply row_spec automatically
apply_row_colors <- function(df, names_caption = NULL,
                             n = 1000, size = 1e3, digits = 5) {
  df <- df %>% 
    as.data.frame() %>% 
    select(where(~ !all(is.na(head(.x, -1))))) %>%
    # Remove na.omit — keep all rows
    (`*`)(size) %>%
    round(digits)
  
  # Apply cell_spec to each cell, handling NAs per column
  df_modified <- df %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(everything(), ~ {
      # Compute colors only on non-NA values, assign na_color to NAs
      colors <- rep("#BBBBBB", length(.))
      non_na_idx <- !is.na(.)
      if (any(non_na_idx)) {
        colors[non_na_idx] <- spec_color2(
          .[non_na_idx],
          option  = "D",
          end     = 0.7,
          palette = paletteer_d("ggsci::amber_material")
        )
      }
      # Format: NA cells show "NA" with gray background, others colored
      cell_spec(
        ifelse(is.na(.), "NA", as.character(.)),
        background = colors,
        color      = "black"
      )
    })) %>%
    t() %>%
    as.data.frame()
  
  # Create a kable table with the modified data frame
  table <- df_modified %>%
    kable(format   = "html", align = "c",
          caption  = names_caption,
          escape   = FALSE) %>%
    kable_styling(full_width = FALSE, c("striped", "condensed"))
  
  return(table)
}
#names_var <- c("catx1","catx2","catx3")