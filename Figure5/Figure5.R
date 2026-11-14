
# ============================================================
# FIGURE 5: Plasma / Proteome ratios by cohort
# ============================================================

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. File location
# ------------------------------------------------------------

folder <- "Your Input Path Here"
data_file <- file.path(folder, "dataFigure5.csv")

if (!file.exists(data_file)) {
  stop("No encuentro el archivo: ", data_file)
}

# ------------------------------------------------------------
# 2. Import data
# ------------------------------------------------------------

dat <- read.csv(
  data_file,
  check.names = FALSE,
  na.strings = c("ND", "")
)

# ------------------------------------------------------------
# 3. Amino acid groups
# ------------------------------------------------------------

non_eaa <- c(
  "Ala", "Asp", "Glu", "Asn", "Pro", "Gln",
  "Arg", "Ser", "Gly", "Cys", "Tyr"
)

eaa <- c(
  "Phe", "His", "Ile", "Leu", "Lys",
  "Met", "Thr", "Trp", "Val"
)

# ------------------------------------------------------------
# 4. Cohort definitions
# ------------------------------------------------------------

cohort_columns <- c(
  "Healthy_Human",
  "DM_Human",
  "MetabolicSyndrome_Human",
  "Dyslipidemia_Human",
  "Hypertension_Human"
)

cohort_names <- c(
  "Healthy_Human" = "Healthy control",
  "DM_Human" = "Diabetes mellitus",
  "MetabolicSyndrome_Human" = "Metabolic syndrome",
  "Dyslipidemia_Human" = "Dyslipidemia",
  "Hypertension_Human" = "Hypertension"
)

cohort_order <- unname(cohort_names)

# ------------------------------------------------------------
# 5. Cohort colors and shapes
# ------------------------------------------------------------

cohort_colors <- c(
  "Healthy control" = "#111111",
  "Diabetes mellitus" = "#1456D9",
  "Metabolic syndrome" = "#159B39",
  "Dyslipidemia" = "#F2A000",
  "Hypertension" = "#E3271A"
)

cohort_shapes <- c(
  "Healthy control" = 16,
  "Diabetes mellitus" = 15,
  "Metabolic syndrome" = 17,
  "Dyslipidemia" = 18,
  "Hypertension" = 16
)

# Horizontal offsets

cohort_offsets <- c(
  "Healthy control" = -0.30,
  "Diabetes mellitus" = -0.15,
  "Metabolic syndrome" = 0,
  "Dyslipidemia" = 0.15,
  "Hypertension" = 0.30
)

# ------------------------------------------------------------
# 6. Check required columns
# ------------------------------------------------------------

required_columns <- c("AA", cohort_columns)

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if (length(missing_columns) > 0) {
  stop(
    "Faltan columnas en dataFigure5.csv: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(dat$AA)) {
  stop("Debe existir una sola fila por aminoácido.")
}

# Check amino acids

all_amino_acids <- c(non_eaa, eaa)

missing_aa <- setdiff(all_amino_acids, dat$AA)

if (length(missing_aa) > 0) {
  stop(
    "Faltan aminoácidos en dataFigure5.csv: ",
    paste(missing_aa, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 7. Convert cohort values to numeric
# ------------------------------------------------------------

missing_count <- 0

for (column in cohort_columns) {
  
  values <- suppressWarnings(
    as.numeric(as.character(dat[[column]]))
  )
  
  missing_count <- missing_count +
    sum(values == 0, na.rm = TRUE)
  
  # Treat zero as missing
  values[values == 0] <- NA_real_
  
  dat[[column]] <- values
}

if (missing_count > 0) {
  message(
    missing_count,
    " valores cero tratados como ausentes."
  )
}

# ------------------------------------------------------------
# 8. Convert data to long format
# ------------------------------------------------------------

plot_data <- do.call(
  rbind,
  lapply(cohort_columns, function(column) {
    
    data.frame(
      AA = dat$AA,
      Cohort = unname(cohort_names[column]),
      Ratio = dat[[column]],
      stringsAsFactors = FALSE
    )
    
  })
)

plot_data$Cohort <- factor(
  plot_data$Cohort,
  levels = cohort_order
)

plot_data$AA <- as.character(plot_data$AA)

# Remove missing and non-positive values

plot_data <- subset(
  plot_data,
  !is.na(Ratio) &
    is.finite(Ratio) &
    Ratio > 0
)

# ------------------------------------------------------------
# 9. Function to create panels
# ------------------------------------------------------------

make_panel <- function(
    amino_acids,
    panel_fill,
    show_y_axis = TRUE
) {
  
  d <- subset(
    plot_data,
    AA %in% amino_acids
  )
  
  d$x <- match(
    d$AA,
    amino_acids
  ) + cohort_offsets[as.character(d$Cohort)]
  
  n <- length(amino_acids)
  
  ggplot() +
    
    # Background
    
    annotate(
      "rect",
      xmin = 0.5,
      xmax = n + 0.5,
      ymin = 0.25,
      ymax = 20,
      fill = panel_fill,
      colour = NA
    ) +
    
    # Reference line at ratio = 1
    
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "#555555",
      linewidth = 0.7
    ) +
    
    # Cohort points
    
    geom_point(
      data = d,
      aes(
        x = x,
        y = Ratio,
        colour = Cohort,
        shape = Cohort
      ),
      size = 3.1
    ) +
    
    # X axis
    
    scale_x_continuous(
      breaks = seq_len(n),
      labels = amino_acids,
      limits = c(0.5, n + 0.5),
      expand = c(0, 0)
    ) +
    
    # Logarithmic Y axis
    
    scale_y_log10(
      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20),
      labels = c(
        "0.25", "0.5", "1", "2",
        "5", "10", "20"
      ),
      limits = c(0.25, 20)
    ) +
    
    # Cohort colors without legend
    
    scale_colour_manual(
      values = cohort_colors,
      guide = "none"
    ) +
    
    # Cohort shapes without legend
    
    scale_shape_manual(
      values = cohort_shapes,
      guide = "none"
    ) +
    
    coord_cartesian(
      ylim = c(0.25, 20)
    ) +
    
    # Axis labels
    
    labs(
      x = NULL,
      y = if (show_y_axis) {
        "RQAA (Plasma) / RQAA (Proteome)"
      } else {
        NULL
      }
    ) +
    
    theme_classic(
      base_size = 15
    ) +
    
    theme(
      
      panel.background = element_blank(),
      panel.grid = element_blank(),
      
      axis.line = element_line(
        colour = "#333333"
      ),
      
      # Amino acid labels: size 15
      
      axis.text.x = element_text(
        size = 14,
        colour = "#111111"
      ),
      
      # Y axis tick labels: size 12
      
      axis.text.y = if (show_y_axis) {
        element_text(
          size = 14,
          colour = "#111111"
        )
      } else {
        element_blank()
      },
      
      axis.ticks.y = if (show_y_axis) {
        element_line()
      } else {
        element_blank()
      },
      
      axis.line.y = if (show_y_axis) {
        element_line(colour = "#333333")
      } else {
        element_blank()
      },
      
      # Y axis title: size 16
      
      axis.title.y = if (show_y_axis) {
        element_text(
          size = 16,
          margin = margin(r = 8)
        )
      } else {
        element_blank()
      },
      
      # Remove all legends
      
      legend.position = "none",
      
      plot.margin = margin(4, 2, 4, 2)
    )
}

# ------------------------------------------------------------
# 10. Create panels
# ------------------------------------------------------------

p_non_eaa <- make_panel(
  amino_acids = non_eaa,
  panel_fill = "#E4F3FC",
  show_y_axis = TRUE
)

p_eaa <- make_panel(
  amino_acids = eaa,
  panel_fill = "#FDEAEA",
  show_y_axis = FALSE
)


# ------------------------------------------------------------
# 11. Combine panels
# ------------------------------------------------------------

figure <- (p_non_eaa + p_eaa) +
  plot_layout(
    widths = c(
      length(non_eaa),
      length(eaa)
    )
  )

# ------------------------------------------------------------
# 12. Export high-resolution TIFF
# ------------------------------------------------------------

ggsave(
  filename = file.path(folder, "Figure5.tiff"),
  plot = figure,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

# ------------------------------------------------------------
# 13. Export TIFF
# ------------------------------------------------------------

ggsave(
  filename = file.path(folder, "Figure5.tiff"),
  plot = figure,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 15. Display figure
# ------------------------------------------------------------

print(figure)

message(
  "Figure 5 saved successfully in: ",
  folder
)
