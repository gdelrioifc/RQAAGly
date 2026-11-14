
# ============================================================
# FIGURE 7: Environmental impacts by RQAA group
# Welch t-tests + Holm correction
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# 1. File location
# ------------------------------------------------------------

folder <- "Your Path Here"

normal_file <- file.path(folder, "Figure7A.tsv")
excess_file <- file.path(folder, "Figure7B.tsv")

if (!file.exists(normal_file) || !file.exists(excess_file)) {
  stop("Put Figure7A.tsv and Figure7B.tsv in: ", folder)
}

# ------------------------------------------------------------
# 2. Read input files
# ------------------------------------------------------------

required <- c(
  "Index", "Description", "GHE", "LandUse", "WaterUse"
)

read_group <- function(path, group_name) {
  
  x <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "ND"),
    quote = "\""
  )
  
  missing_cols <- setdiff(required, names(x))
  
  if (length(missing_cols)) {
    stop(
      basename(path),
      " is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  for (nm in c("GHE", "LandUse", "WaterUse")) {
    
    x[[nm]] <- suppressWarnings(
      as.numeric(x[[nm]])
    )
    
    if (any(!is.finite(x[[nm]]) | x[[nm]] <= 0)) {
      stop(
        basename(path),
        " contains invalid values in ",
        nm,
        "; log-scale plotting requires positive measurements."
      )
    }
  }
  
  data.frame(
    Index = x$Index,
    Description = x$Description,
    GHE = x$GHE,
    LandUse = x$LandUse,
    WaterUse = x$WaterUse,
    Group = group_name,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 3. Load groups
# ------------------------------------------------------------

normal <- read_group(
  normal_file,
  "RQAA <= 2.0"
)

excess <- read_group(
  excess_file,
  "RQAA > 2.0"
)

if (nrow(normal) != 412L || nrow(excess) != 128L) {
  warning(
    "The TSV row counts differ from the expected 412 normal and 128 excess foods."
  )
}

foods <- rbind(normal, excess)

if (anyDuplicated(foods$Index)) {
  warning(
    "Some food Index values occur more than once across the two TSV files."
  )
}

foods$Group <- factor(
  foods$Group,
  levels = c("RQAA <= 2.0", "RQAA > 2.0")
)

# ------------------------------------------------------------
# 4. Environmental metrics
# ------------------------------------------------------------

metrics <- list(
  
  GHE = list(
    title = "Greenhouse gas emissions"
  ),
  
  LandUse = list(
    title = "Land use"
  ),
  
  WaterUse = list(
    title = "Water use"
  )
  
)

# ------------------------------------------------------------
# 5. Long format
# ------------------------------------------------------------

long <- do.call(
  rbind,
  lapply(names(metrics), function(nm) {
    
    data.frame(
      Index = foods$Index,
      Group = foods$Group,
      Metric = nm,
      Value = foods[[nm]],
      stringsAsFactors = FALSE
    )
    
  })
)

# ------------------------------------------------------------
# 6. Statistical tests and Holm correction
# ------------------------------------------------------------

p_values <- setNames(
  numeric(length(metrics)),
  names(metrics)
)

# Calculate all Welch tests before plotting

for (nm in names(metrics)) {
  
  d <- long[
    long$Metric == nm,
    ,
    drop = FALSE
  ]
  
  test <- t.test(
    log10(Value) ~ Group,
    data = d
  )
  
  p_values[nm] <- test$p.value
}

# Apply Holm correction across the three comparisons

p_adjusted <- p.adjust(
  p_values,
  method = "holm"
)

# Display statistical results

statistical_results <- data.frame(
  Metric = names(p_values),
  p_Welch = as.numeric(p_values),
  p_Holm = as.numeric(p_adjusted),
  Significant_Holm = p_adjusted < 0.05
)

print(statistical_results)

summary_rows <- list()

# ------------------------------------------------------------
# 7. Create individual plots
# ------------------------------------------------------------

plots <- lapply(names(metrics), function(nm) {
  
  d <- long[
    long$Metric == nm,
    ,
    drop = FALSE
  ]
  
  # Holm-adjusted p-value with subscript
  
  p_label <- if (p_adjusted[nm] < 0.001) {
    
    "italic(p)[adj] < 0.001"
    
  } else {
    
    paste0(
      "italic(p)[adj] == ",
      formatC(
        p_adjusted[nm],
        digits = 3,
        format = "f"
      )
    )
  }
  
  # Descriptive statistics
  
  group_summary <- do.call(
    rbind,
    lapply(levels(foods$Group), function(g) {
      
      v <- d$Value[d$Group == g]
      
      data.frame(
        Metric = nm,
        Group = g,
        n = length(v),
        median = median(v),
        Q1 = unname(quantile(v, 0.25)),
        Q3 = unname(quantile(v, 0.75))
      )
      
    })
  )
  
  summary_rows[[nm]] <<- group_summary
  
  # Position of statistical annotation
  
  max_value <- max(d$Value)
  
  bracket_y <- max_value * 2.0
  text_y <- max_value * 2.9
  
  # Plot
  
  ggplot(
    d,
    aes(
      x = Group,
      y = Value,
      color = Group,
      shape = Group
    )
  ) +
    
    # Individual food observations
    
    geom_jitter(
      width = 0.18,
      height = 0,
      alpha = 0.55,
      size = 1.15,
      show.legend = FALSE
    ) +
    
    # Boxplots
    
    geom_boxplot(
      aes(group = Group),
      width = 0.28,
      outlier.shape = NA,
      fill = NA,
      color = "#303030",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    
    # Statistical bracket
    
    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y,
      linewidth = 0.35,
      color = "#303030"
    ) +
    
    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = bracket_y * 0.92,
      yend = bracket_y,
      linewidth = 0.35,
      color = "#303030"
    ) +
    
    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = bracket_y * 0.92,
      yend = bracket_y,
      linewidth = 0.35,
      color = "#303030"
    ) +
    
    # Holm-adjusted p-value with subscript
    
    annotate(
      "text",
      x = 1.5,
      y = text_y,
      label = p_label,
      parse = TRUE,
      size = 3.5,
      color = "#303030"
    ) +
    
    # Colors without legend
    
    scale_color_manual(
      values = c(
        "RQAA <= 2.0" = "#3F8CFF",
        "RQAA > 2.0" = "#222222"
      ),
      guide = "none"
    ) +
    
    # Shapes without legend
    
    scale_shape_manual(
      values = c(
        "RQAA <= 2.0" = 16,
        "RQAA > 2.0" = 17
      ),
      guide = "none"
    ) +
    
    # X-axis labels
    
    # X-axis labels without sample sizes
    
    scale_x_discrete(
      labels = c(
        "RQAA <= 2.0" = "RQAA <= 2.0",
        "RQAA > 2.0" = "RQAA > 2.0"
      )
    ) +
    
    # Logarithmic Y axis
    
    scale_y_log10() +
    
    coord_cartesian(
      ylim = c(
        min(d$Value) * 0.65,
        max_value * 3.5
      ),
      clip = "off"
    ) +
    
    # Titles: only environmental metric names
    
    labs(
      title = metrics[[nm]]$title,
      x = NULL,
      y = NULL
    ) +
    
    theme_classic(
      base_size = 15
    ) +
    
    theme(
      
      # Panel titles
      
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      
      # X-axis group labels
      
      axis.text.x = element_text(
        size = 15,
        colour = "#222222"
      ),
      
      # Y-axis numerical labels
      
      axis.text.y = element_text(
        size = 15,
        colour = "#222222"
      ),
      
      # No axis titles
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      # No legends
      
      legend.position = "none",
      
      plot.margin = margin(
        8, 8, 8, 8
      )
    )
  
})

# ------------------------------------------------------------
# 8. Combine three panels
# ------------------------------------------------------------

p <- wrap_plots(
  plots,
  nrow = 1
)

# ------------------------------------------------------------
# 9. Export summary and statistical tests
# ------------------------------------------------------------

summary_table <- do.call(
  rbind,
  summary_rows
)

# Original Welch p-values

summary_table$p_value_log10_welch <- unname(
  p_values[summary_table$Metric]
)

# Holm-adjusted p-values

summary_table$p_value_Holm <- unname(
  p_adjusted[summary_table$Metric]
)

# Significance after Holm correction

summary_table$Significant_Holm <- (
  summary_table$p_value_Holm < 0.05
)

write.csv(
  summary_table,
  file.path(
    folder,
    "figure7_summary_and_tests.csv"
  ),
  row.names = FALSE
)

# Export statistical results

write.csv(
  statistical_results,
  file.path(
    folder,
    "figure7_statistical_tests_Holm.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Data verification
# ------------------------------------------------------------

n_identical <- sum(
  foods$GHE == foods$WaterUse
)

message(
  "Foods: ",
  nrow(normal),
  " with RQAA <= 2.0; ",
  nrow(excess),
  " with RQAA > 2.0."
)

message(
  "Welch p-values: ",
  paste(
    names(p_values),
    format(p_values, digits = 4),
    collapse = "; "
  )
)

message(
  "Holm-adjusted p-values: ",
  paste(
    names(p_adjusted),
    format(p_adjusted, digits = 4),
    collapse = "; "
  )
)

message(
  "GHE and WaterUse are numerically identical for ",
  n_identical,
  " of ",
  nrow(foods),
  " foods. Verify source data if unexpected."
)

# ------------------------------------------------------------
# 11. Export TIFF
# ------------------------------------------------------------

ggsave(
  file.path(
    folder,
    "figure7_environmental_comparison.tiff"
  ),
  p,
  width = 14,
  height = 5.5,
  units = "in",
  dpi = 320,
  bg = "white",
  compression = "lzw"
)

# ------------------------------------------------------------
# 12. Export PNG
# ------------------------------------------------------------

ggsave(
  file.path(
    folder,
    "figure7_environmental_comparison.png"
  ),
  p,
  width = 14,
  height = 5.5,
  units = "in",
  dpi = 320,
  bg = "white"
)

# ------------------------------------------------------------
# 13. Export PDF
# ------------------------------------------------------------

ggsave(
  file.path(
    folder,
    "figure7_environmental_comparison.pdf"
  ),
  p,
  width = 14,
  height = 5.5,
  units = "in",
  bg = "white"
)

# ------------------------------------------------------------
# 14. Display
# ------------------------------------------------------------

print(p)

message(
  "Figure 7 saved successfully in: ",
  folder
)
