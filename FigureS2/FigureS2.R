
# ============================================================
# FIGURE 5: Heatmap - Disease cohorts relative to healthy control
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. File locations
# ------------------------------------------------------------

data_folder <- "Your Input Path Here"

output_folder <- "Your Output Path Here"

stats_file <- file.path(
  data_folder,
  "stastistics_dataFigure5.csv"
)

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

if (!file.exists(stats_file)) {
  stop("Cannot find the statistics CSV: ", stats_file)
}

# ------------------------------------------------------------
# 2. Load statistics
# ------------------------------------------------------------

stats <- read.csv(
  stats_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Amino Acid",
  "Ratio RQAA(free cells)/RQAA(proteome)",
  "p-value",
  "Condition"
)

missing_columns <- setdiff(
  required_columns,
  names(stats)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. Amino acid order
# ------------------------------------------------------------

non_eaa <- c(
  "Ala", "Asp", "Glu", "Asn", "Pro", "Gln",
  "Arg", "Ser", "Gly", "Cys", "Tyr"
)

eaa <- c(
  "Phe", "His", "Ile", "Leu", "Lys",
  "Met", "Thr", "Trp", "Val"
)

aa_order <- c(non_eaa, eaa)

# ------------------------------------------------------------
# 4. Disease cohort order
# ------------------------------------------------------------

disease_order <- c(
  "Dementia",
  "Metabolic syndrome",
  "Dyslipidemia",
  "Hypertension"
)

# ------------------------------------------------------------
# 5. Prepare data
# ------------------------------------------------------------

heat_data <- data.frame(
  AA = stats[["Amino Acid"]],
  Condition = stats[["Condition"]],
  Ratio = suppressWarnings(
    as.numeric(
      stats[["Ratio RQAA(free cells)/RQAA(proteome)"]]
    )
  ),
  p_value = suppressWarnings(
    as.numeric(stats[["p-value"]])
  ),
  stringsAsFactors = FALSE
)

heat_data$Disease <- sub(
  "^Healthy/",
  "",
  heat_data$Condition
)

heat_data$Disease[
  heat_data$Disease == "DM"
] <- "Dementia"

heat_data$Disease[
  heat_data$Disease == "Metabolic Syndrome"
] <- "Metabolic syndrome"

heat_data <- heat_data[
  heat_data$AA %in% aa_order &
    heat_data$Disease %in% disease_order,
]

if (any(heat_data$Ratio <= 0, na.rm = TRUE)) {
  warning(
    "Non-positive ratios treated as missing."
  )
  
  heat_data$Ratio[
    heat_data$Ratio <= 0
  ] <- NA_real_
}

# ------------------------------------------------------------
# 6. Complete amino acid x disease grid
# ------------------------------------------------------------

all_cells <- expand.grid(
  AA = aa_order,
  Disease = disease_order,
  stringsAsFactors = FALSE
)

heat_data <- merge(
  all_cells,
  heat_data[, c(
    "AA", "Disease", "Ratio", "p_value"
  )],
  by = c("AA", "Disease"),
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# 7. Amino acid classes
# ------------------------------------------------------------

heat_data$Class <- ifelse(
  heat_data$AA %in% non_eaa,
  "Non-essential amino acids",
  "Essential amino acids"
)

heat_data$Class <- factor(
  heat_data$Class,
  levels = c(
    "Non-essential amino acids",
    "Essential amino acids"
  )
)

# ------------------------------------------------------------
# 8. Factor ordering
# ------------------------------------------------------------

heat_data$AA <- factor(
  heat_data$AA,
  levels = rev(aa_order)
)

heat_data$Disease <- factor(
  heat_data$Disease,
  levels = disease_order
)

# ------------------------------------------------------------
# 9. Statistical significance
# ------------------------------------------------------------

heat_data$Significance <- ifelse(
  is.na(heat_data$p_value),
  "",
  ifelse(
    heat_data$p_value < 0.001,
    "***",
    ifelse(
      heat_data$p_value < 0.01,
      "**",
      ifelse(
        heat_data$p_value < 0.05,
        "*",
        "ns"
      )
    )
  )
)

# ------------------------------------------------------------
# 10. Cell labels
# ------------------------------------------------------------

heat_data$CellLabel <- ifelse(
  is.na(heat_data$Ratio),
  "",
  sprintf("%.2f", heat_data$Ratio)
)

has_stars <- (
  !is.na(heat_data$Ratio) &
    heat_data$Significance != ""
)

heat_data$CellLabel[has_stars] <- paste0(
  heat_data$CellLabel[has_stars],
  "\n",
  heat_data$Significance[has_stars]
)

# ------------------------------------------------------------
# 11. Log2 ratio
# ------------------------------------------------------------

heat_data$Log2Ratio <- log2(
  heat_data$Ratio
)

max_abs <- max(
  abs(heat_data$Log2Ratio),
  na.rm = TRUE
)

if (!is.finite(max_abs) || max_abs == 0) {
  max_abs <- 0.1
}

# ------------------------------------------------------------
# 12. Text contrast
# ------------------------------------------------------------

heat_data$TextColor <- ifelse(
  !is.na(heat_data$Log2Ratio) &
    abs(heat_data$Log2Ratio) >= 0.35,
  "white",
  "#222222"
)

# ------------------------------------------------------------
# 13. Create heatmap
# ------------------------------------------------------------

heatmap <- ggplot(
  heat_data,
  aes(
    x = Disease,
    y = AA,
    fill = Log2Ratio
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +
  
  geom_text(
    aes(
      label = CellLabel,
      colour = TextColor
    ),
    size = 3.1,
    lineheight = 0.9
  ) +
  
  facet_grid(
    Class ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_colour_identity() +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    na.value = "#E6E6E6",
    guide = "none"
  ) +
  
  # ----------------------------------------------------------
# 14. Theme and formatting
# ----------------------------------------------------------

theme_minimal(base_size = 15) +
  
  theme(
    
    panel.grid = element_blank(),
    
    # Remove axis titles AA and Disease
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Disease names: size 15
    
    axis.text.x = element_text(
      size = 15,
      colour = "#333333",
      angle = 0,
      hjust = 0.5
    ),
    
    # Amino acid names: size 15
    
    axis.text.y = element_text(
      size = 15,
      colour = "#333333"
    ),
    
    # Facet labels
    
    strip.placement = "outside",
    
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "#F0F0F0",
      colour = NA
    ),
    
    # Panel spacing
    
    panel.spacing.y = grid::unit(
      0.4,
      "lines"
    ),
    
    # Remove legends completely
    
    legend.position = "none",
    
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

# ------------------------------------------------------------
# 15. Export TIFF
# ------------------------------------------------------------

ggsave(
  file.path(
    output_folder,
    "figure5.1_heatmap.tiff"
  ),
  plot = heatmap,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 17. Display
# ------------------------------------------------------------

print(heatmap)

message(
  "PNG and PDF saved successfully in: ",
  output_folder
)

