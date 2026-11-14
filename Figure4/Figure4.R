# ============================================================
# FIGURE 4 - FREE-POOL-TO-PROTEOME AMINO ACID RATIOS
# Group medians of species-level medians + interquartile range
# ============================================================

library(tidyverse)

# 1. File paths
input_file <- "dataFigure4.csv"

output_dir <- "Your OutPut Path Here"

# 2. Import summary CSV
summary <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 3. Prepare plotting variables
summary <- summary %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Plants", "Bacteria", "Fungi", "Animals")
    ),
    Class = factor(
      Class,
      levels = c("Non-EAA", "EAA")
    ),
    median_label = sprintf("%.2f", median_between_species),
    n_label = paste0("n = ", n_species),
    
    # Explicit x coordinates for labels.
    # Bars are dodged by +/- 0.195 around each group center.
    # Add +0.10 to shift BOTH median labels slightly to the right.
    x_group = as.numeric(Group),
    x_median_label = x_group +
      if_else(Class == "Non-EAA", -0.195, 0.195) + 0.10,
    
    median_label_y = if_else(
      median_log10 >= 0,
      median_log10 + 0.085,
      median_log10 - 0.085
    ),
    n_label_y = Q3_log10 + 0.25
  )

# 4. Plot
# Bar height: log10 of the group median of species-level medians.
# Error bars: Q1-Q3 across species.
# Dashed line at 0: ratio = 1.

dodge <- position_dodge(width = 0.78)

p <- ggplot(
  summary,
  aes(
    x = Group,
    y = median_log10,
    fill = Class,
    group = Class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.6
  ) +
  geom_col(
    position = dodge,
    width = 0.64,
    color = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = Q1_log10,
      ymax = Q3_log10
    ),
    position = dodge,
    width = 0.16,
    linewidth = 0.55
  ) +
  
  # Non-EAA = blue; EAA = red. No Class legend.
  scale_fill_manual(
    values = c(
      "Non-EAA" = "#0066A8",
      "EAA" = "#C62828"
    )
  ) +
  # Median labels use explicit x coordinates, so the rightward shift
  # is applied directly rather than relying on nudge_x with dodge.
  geom_text(
    aes(
      x = x_median_label,
      y = median_label_y,
      label = median_label
    ),
    size = 4.5,
    fontface = "bold",
    color = "black",
    vjust = 0.5,
    inherit.aes = FALSE,
    data = summary
  ) +
  # Species counts remain centered above each IQR.
  geom_text(
    aes(
      y = n_label_y,
      label = n_label
    ),
    position = dodge,
    size = 3.6,
    color = "gray30"
  ) +
  scale_x_discrete(
    limits = levels(summary$Group),
    expand = expansion(add = 0.55)
  ) +
  scale_y_continuous(
    breaks = seq(-0.5, 2.0, by = 0.5),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.18, 0.24))
  ) +
  
  labs(
    x = NULL,
    y = expression(log[10] * " (Free-pool RQAA / Proteome RQAA)")
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "gray30"
    ),
    axis.title.y = element_text(
      size = 15,
      face = "bold",
      margin = margin(r = 15)
    ),
    axis.text.y = element_text(
      size = 15,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.caption = element_text(
      size = 8.5,
      color = "gray35",
      hjust = 0,
      lineheight = 1.15
    ),
    plot.caption.position = "plot",
    plot.margin = margin(20, 18, 14, 25)
  )

print(p)

# 5. Export high-resolution TIFF
ggsave(
  filename = file.path(output_dir, "Figure4_species_medians_IQR_final.tiff"),
  plot = p,
  width = 10,
  height = 6.5,
  dpi = 600,
  bg = "white"
)


# 6. Print summary values
summary %>%
  select(
    Group,
    Class,
    n_species,
    median_between_species,
    Q1,
    Q3,
    IQR
  ) %>%
  print()
