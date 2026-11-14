
# ============================================================
# FIGURE 6: Food values below and above RQAA = 2
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. File location
# ------------------------------------------------------------

folder <- "Your Path Here"

data_file <- file.path(
  folder,
  "dataFigure6_WithFoodNames.csv"
)

if (!file.exists(data_file)) {
  stop("CSV not found: ", data_file)
}

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

dat <- read.csv(
  data_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "ND"),
  fileEncoding = "UTF-8-BOM"
)

if (!all(c("Index", "AA") %in% names(dat)) ||
    ncol(dat) <= 2L) {
  stop("Expected Index, AA, and food value columns in the CSV.")
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

if (!all(aa_order %in% dat$AA)) {
  stop("The CSV must contain all 20 expected amino acids.")
}

# ------------------------------------------------------------
# 4. Prepare values
# ------------------------------------------------------------

value_matrix <- as.matrix(
  dat[, 3:ncol(dat), drop = FALSE]
)

suppressWarnings(
  storage.mode(value_matrix) <- "numeric"
)

long <- data.frame(
  AA = rep(
    as.character(dat$AA),
    each = ncol(value_matrix)
  ),
  RQAA = as.vector(t(value_matrix)),
  stringsAsFactors = FALSE
)

# Exclude missing, zero and non-positive values

long$Observed <- (
  is.finite(long$RQAA) &
    long$RQAA > 0
)

long$Class <- ifelse(
  long$AA %in% eaa,
  "Essential",
  "Non-essential"
)

n_foods <- ncol(value_matrix)

# ------------------------------------------------------------
# 5. Calculate percentages using RQAA = 2
# ------------------------------------------------------------

summary_data <- do.call(
  rbind,
  lapply(aa_order, function(aa) {
    
    x <- long$RQAA[
      long$AA == aa & long$Observed
    ]
    
    # Values below, above and exactly equal to 2
    
    n_below <- sum(x < 2)
    n_above <- sum(x > 2)
    n_equal <- sum(x == 2)
    
    # Exact ties are excluded from the denominator
    
    n_compared <- n_below + n_above
    
    data.frame(
      AA = aa,
      
      Class = if (aa %in% eaa) {
        "Essential"
      } else {
        "Non-essential"
      },
      
      n_observed = length(x),
      n_ND = n_foods - length(x),
      
      n_below_2 = n_below,
      
      percent_below_2 = if (n_compared > 0) {
        100 * n_below / n_compared
      } else {
        NA_real_
      },
      
      n_equal_2 = n_equal,
      
      n_above_2 = n_above,
      
      percent_above_2 = if (n_compared > 0) {
        100 * n_above / n_compared
      } else {
        NA_real_
      }
    )
  })
)

summary_data$AA <- factor(
  summary_data$AA,
  levels = rev(aa_order)
)

summary_data$Class <- factor(
  summary_data$Class,
  levels = c("Non-essential", "Essential")
)

# ------------------------------------------------------------
# 6. Export summary CSV for RQAA = 2
# ------------------------------------------------------------

write.csv(
  summary_data,
  file.path(
    folder,
    "figure6_1_summary_below_above_RQAA2.csv"
  ),
  row.names = FALSE,
  na = "ND"
)

# ------------------------------------------------------------
# 7. Prepare diverging bars
# ------------------------------------------------------------

bar_data <- rbind(
  
  data.frame(
    AA = summary_data$AA,
    Class = summary_data$Class,
    Side = "Below 2",
    Percent = -summary_data$percent_below_2
  ),
  
  data.frame(
    AA = summary_data$AA,
    Class = summary_data$Class,
    Side = "Above 2",
    Percent = summary_data$percent_above_2
  )
)

# Percentage labels

bar_data$label <- ifelse(
  is.na(bar_data$Percent),
  "",
  paste0(
    round(abs(bar_data$Percent)),
    "%"
  )
)

# Label positions:
# Below 2: outside the left bar
# Above 2: outside the right bar
# Zero values: slightly to the right of the center line

bar_data$label_x <- ifelse(
  bar_data$Percent < 0,
  bar_data$Percent - 3,
  ifelse(
    bar_data$Percent > 0,
    bar_data$Percent + 3,
    3
  )
)

bar_data$hjust <- ifelse(
  bar_data$Percent < 0,
  1,
  0
)

# ------------------------------------------------------------
# 8. Handle amino acids without observations
# ------------------------------------------------------------

nd_aa <- summary_data$AA[
  summary_data$n_observed == 0
]

nd_data <- data.frame(
  AA = factor(
    nd_aa,
    levels = rev(aa_order)
  ),
  x = 0,
  label = "ND"
)

# ------------------------------------------------------------
# 9. Colors
# ------------------------------------------------------------

colors <- c(
  "Non-essential" = "#2878B5",
  "Essential" = "#D94841"
)

# ------------------------------------------------------------
# 10. Create plot
# ------------------------------------------------------------

p <- ggplot(
  bar_data,
  aes(
    x = Percent,
    y = AA,
    fill = Class
  )
) +
  
  # Bars
  
  geom_col(
    width = 0.68,
    na.rm = TRUE
  ) +
  
  # Central reference line
  
  geom_vline(
    xintercept = 0,
    color = "#303030",
    linewidth = 0.55
  ) +
  
  # Separation between amino acid groups
  
  geom_hline(
    yintercept = length(eaa) + 0.5,
    color = "#999999",
    linewidth = 0.45
  ) +
  
  # Percentage labels
  
  geom_text(
    data = subset(bar_data, label != ""),
    aes(
      x = label_x,
      y = AA,
      label = label,
      hjust = hjust
    ),
    inherit.aes = FALSE,
    size = 4,
    color = "#252525"
  ) +
  
  # ND labels
  
  geom_text(
    data = nd_data,
    aes(
      x = x,
      y = AA,
      label = label
    ),
    inherit.aes = FALSE,
    color = "#666666",
    fontface = "italic",
    size = 3.5
  ) +
  
  # Fill colors
  
  scale_fill_manual(
    values = colors,
    guide = "none"
  ) +
  
  # X axis
  
  scale_x_continuous(
    limits = c(-112, 112),
    breaks = seq(-100, 100, 25),
    labels = function(x) {
      paste0(abs(x), "%")
    },
    expand = expansion(mult = c(0, 0))
  ) +
  
  # Y axis
  
  scale_y_discrete(
    labels = function(x) as.character(x),
    drop = FALSE
  ) +
  
  # Remove axis titles
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Amino acid labels
    
    axis.text.y = element_text(
      face = "plain",
      size = 12,
      colour = "#222222"
    ),
    
    # X-axis percentage labels
    
    axis.text.x = element_text(
      size = 12,
      colour = "#222222"
    ),
    
    # No legend
    
    legend.position = "none",
    
    plot.margin = margin(
      10, 20, 8, 10
    )
  )

# ------------------------------------------------------------
# 11. Display plot
# ------------------------------------------------------------

print(p)

# ------------------------------------------------------------
# 12. Export TIFF
# ------------------------------------------------------------

ggsave(
  file.path(folder, "Figure6_RQAA2.tiff"),
  plot = p,
  width = 12,
  height = 9,
  units = "in",
  dpi = 320,
  bg = "white",
  compression = "lzw"
)


# ------------------------------------------------------------
# 15. Print overall percentages
# ------------------------------------------------------------

overall_results <- do.call(
  rbind,
  lapply(
    c("Non-essential", "Essential"),
    function(cl) {
      
      x <- summary_data[
        summary_data$Class == cl,
      ]
      
      n_below <- sum(x$n_below_2)
      n_above <- sum(x$n_above_2)
      
      data.frame(
        Class = cl,
        n_below_2 = n_below,
        n_above_2 = n_above,
        percent_below_2 = 100 * n_below /
          (n_below + n_above),
        percent_above_2 = 100 * n_above /
          (n_below + n_above)
      )
    }
  )
)

print(overall_results)

message(
  "Figure 6 (RQAA = 2) and summary CSV saved successfully in: ",
  folder
)

