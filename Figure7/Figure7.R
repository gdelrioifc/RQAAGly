# ============================================================
# FIGURE 7: Environmental impacts by RQAA group
# Two-sided Welch's t-tests on log10-transformed measurements
# No multiple-comparison correction
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# 1. File location

folder <- "C:/Users/sofii/OneDrive/Documentos/Figura1_RQAA/Figura 7"

normal_file <- file.path(folder, "Figure7A.tsv")
excess_file <- file.path(folder, "Figure7B.tsv")

if (!file.exists(normal_file) || !file.exists(excess_file)) {
  stop(
    "Could not find the input TSV files. Check folder and filenames:\n",
    normal_file, "\n", excess_file
  )
}

# ------------------------------------------------------------
# 2. Read input files
# ------------------------------------------------------------

required <- c("Index", "Description", "GHE", "LandUse", "WaterUse")

read_group <- function(path, group_name) {
  x <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "ND"),
    quote = "\""
  )

  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop(
      basename(path), " is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  for (nm in c("GHE", "LandUse", "WaterUse")) {
    x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))

    # Welch tests are performed on log10 values, so measurements
    # must be finite and strictly positive.
    if (any(!is.finite(x[[nm]]) | x[[nm]] <= 0)) {
      stop(
        basename(path), " contains missing, non-numeric, or non-positive ",
        "values in ", nm, ". Check the source TSV before analysis."
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

normal <- read_group(normal_file, "RQAA <= 2.0")
excess <- read_group(excess_file, "RQAA > 2.0")

if (nrow(normal) != 412L || nrow(excess) != 128L) {
  warning(
    "The TSV row counts differ from the expected 412 foods with RQAA <= 2.0 ",
    "and 128 foods with RQAA > 2.0."
  )
}

foods <- rbind(normal, excess)

if (anyDuplicated(foods$Index)) {
  warning("Some food Index values occur more than once across the two TSV files.")
}

foods$Group <- factor(
  foods$Group,
  levels = c("RQAA <= 2.0", "RQAA > 2.0")
)

# ------------------------------------------------------------
# 4. Environmental metrics
# ------------------------------------------------------------

metrics <- list(
  GHE = list(title = "Greenhouse gas emissions"),
  LandUse = list(title = "Land use"),
  WaterUse = list(title = "Water use")
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
# 6. Two-sided Welch's t-tests on log10-transformed values
# ------------------------------------------------------------

# t.test() uses a two-sided alternative by default and Welch's
# unequal-variance correction when var.equal = FALSE.
# No Holm-Bonferroni or other multiple-testing correction is applied.

test_results <- lapply(names(metrics), function(nm) {
  d <- long[long$Metric == nm, , drop = FALSE]

  test <- t.test(
    log10(Value) ~ Group,
    data = d,
    alternative = "two.sided",
    var.equal = FALSE
  )

  data.frame(
    Metric = nm,
    n_RQAA_le_2 = sum(d$Group == "RQAA <= 2.0"),
    n_RQAA_gt_2 = sum(d$Group == "RQAA > 2.0"),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value,
    conf_low_log10 = unname(test$conf.int[1]),
    conf_high_log10 = unname(test$conf.int[2]),
    stringsAsFactors = FALSE
  )
})

statistical_results <- do.call(rbind, test_results)
print(statistical_results)

p_values <- setNames(
  statistical_results$p_value,
  statistical_results$Metric
)

# ------------------------------------------------------------
# 7. Create individual plots
# ------------------------------------------------------------

summary_rows <- list()

format_p_label <- function(p) {
  if (p < 0.001) {
    "italic(p) < 0.001"
  } else {
    paste0(
      "italic(p) == ",
      formatC(p, digits = 3, format = "f")
    )
  }
}

plots <- lapply(names(metrics), function(nm) {
  d <- long[long$Metric == nm, , drop = FALSE]

  p_label <- format_p_label(p_values[nm])

  # Descriptive statistics on the original measurement scale.
  group_summary <- do.call(
    rbind,
    lapply(levels(foods$Group), function(g) {
      v <- d$Value[d$Group == g]

      data.frame(
        Metric = nm,
        Group = g,
        n = length(v),
        mean = mean(v),
        SD = sd(v),
        median = median(v),
        Q1 = unname(quantile(v, 0.25)),
        Q3 = unname(quantile(v, 0.75))
      )
    })
  )

  group_summary$p_value_two_sided_welch_log10 <- p_values[nm]
  summary_rows[[nm]] <<- group_summary

  max_value <- max(d$Value)
  bracket_y <- max_value * 1.65
  text_y <- max_value * 2.15

  ggplot(
    d,
    aes(
      x = Group,
      y = Value,
      color = Group,
      shape = Group
    )
  ) +
    geom_jitter(
      width = 0.18,
      height = 0,
      alpha = 0.55,
      size = 1.15,
      show.legend = FALSE
    ) +
    geom_boxplot(
      aes(group = Group),
      width = 0.28,
      outlier.shape = NA,
      fill = NA,
      color = "#303030",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    # Comparison bracket
    annotate(
      "segment",
      x = 1, xend = 2,
      y = bracket_y, yend = bracket_y,
      linewidth = 0.35,
      color = "#303030"
    ) +
    annotate(
      "segment",
      x = 1, xend = 1,
      y = bracket_y * 0.94, yend = bracket_y,
      linewidth = 0.35,
      color = "#303030"
    ) +
    annotate(
      "segment",
      x = 2, xend = 2,
      y = bracket_y * 0.94, yend = bracket_y,
      linewidth = 0.35,
      color = "#303030"
    ) +
    # Unadjusted p-value from the Welch test for this metric only
    annotate(
      "text",
      x = 1.5,
      y = text_y,
      label = p_label,
      parse = TRUE,
      size = 3.5,
      color = "#303030"
    ) +
    scale_color_manual(
      values = c(
        "RQAA <= 2.0" = "#3F8CFF",
        "RQAA > 2.0" = "#222222"
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "RQAA <= 2.0" = 16,
        "RQAA > 2.0" = 17
      ),
      guide = "none"
    ) +
    scale_x_discrete(
      labels = c(
        "RQAA <= 2.0" = "RQAA <= 2.0",
        "RQAA > 2.0" = "RQAA > 2.0"
      )
    ) +
    scale_y_log10() +
    coord_cartesian(
      ylim = c(min(d$Value) * 0.65, max_value * 2.7),
      clip = "off"
    ) +
    labs(
      title = metrics[[nm]]$title,
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      axis.text.x = element_text(
        size = 15,
        colour = "#222222"
      ),
      axis.text.y = element_text(
        size = 15,
        colour = "#222222"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8)
    )
})

# ------------------------------------------------------------
# 8. Combine the three independent environmental comparisons
# ------------------------------------------------------------

p <- wrap_plots(plots, nrow = 1)

# ------------------------------------------------------------
# 9. Export descriptive summary and unadjusted test results
# ------------------------------------------------------------

summary_table <- do.call(rbind, summary_rows)

write.csv(
  summary_table,
  file.path(folder, "figure7_summary_and_unadjusted_tests.csv"),
  row.names = FALSE
)

write.csv(
  statistical_results,
  file.path(folder, "figure7_welch_t_tests_unadjusted.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Data verification
# ------------------------------------------------------------

n_identical_ghe_water <- sum(foods$GHE == foods$WaterUse)

message(
  "Foods loaded: ", nrow(normal), " with RQAA <= 2.0; ",
  nrow(excess), " with RQAA > 2.0."
)

message("Two-sided Welch t-tests on log10-transformed measurements (unadjusted):")
print(statistical_results)

if (n_identical_ghe_water > 0) {
  message(
    "Data check: GHE and WaterUse are numerically identical for ",
    n_identical_ghe_water, " of ", nrow(foods),
    " foods. This is a source-data check only; the metrics are tested ",
    "separately and are not compared against each other."
  )
}

# ------------------------------------------------------------
# 11. Export Figure 7
# ------------------------------------------------------------

ggsave(
  file.path(folder, "Figure7.tiff"),
  p,
  width = 14,
  height = 5.5,
  units = "in",
  dpi = 320,
  bg = "white",
  compression = "lzw"
)


# Display figure in RStudio
print(p)

message("Figure 7 saved successfully in: ", folder)
