
# ============================================================
# Relative amino-acid quantities across taxonomic groups
# ============================================================

# Install ggplot2 once, if needed:
# install.packages("ggplot2")

library(ggplot2)

# 1. Read the summary CSV
data <- read.csv(
  "C:dataFigureS1.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. Rename columns
plot_data <- data[, c(
  "Taxonomic group",
  "Amino-acid class",
  "Valid measurements (n)",
  "RQAA > 1 (%)",
  "RQAA < 1 or = 1 (%)"
)]

names(plot_data) <- c(
  "Group", "Class", "n", "Above", "At_or_below"
)

# 3. Set the group and class order
plot_data$Group <- factor(
  plot_data$Group,
  levels = c(
    "Animals", "Plants", "Fungi",
    "Bacteria", "All taxonomic groups"
  ),
  labels = c(
    "Animals", "Plants", "Fungi",
    "Bacteria", "All groups"
  )
)

plot_data$Class <- factor(
  plot_data$Class,
  levels = c("EAA", "non-EAA")
)

# 4. Create one row per bar segment
plot_long <- rbind(
  data.frame(
    Group = plot_data$Group,
    Class = plot_data$Class,
    n = plot_data$n,
    Category = "RQAA > 1",
    Percentage = plot_data$Above
  ),
  data.frame(
    Group = plot_data$Group,
    Class = plot_data$Class,
    n = plot_data$n,
    Category = "RQAA \u2264 1",
    Percentage = plot_data$At_or_below
  )
)

plot_long$Category <- factor(
  plot_long$Category,
  levels = c("RQAA > 1", "RQAA \u2264 1")
)

# 5. Create positions for paired bars and gaps between groups
group_order <- c(
  "Animals", "Plants", "Fungi", "Bacteria", "All groups"
)

group_centers <- c(1, 2.5, 4, 5.5, 7.5)

plot_long$Group <- factor(
  plot_long$Group,
  levels = group_order
)

plot_long$X <- rep(
  group_centers[match(as.character(plot_long$Group), group_order)],
  each = 1
)

# EAA left; non-EAA right within each pair
plot_long$X <- plot_long$X +
  ifelse(plot_long$Class == "EAA", -0.23, 0.23)

# 6. Percentage labels
plot_long$Label <- paste0(
  sprintf("%.1f", plot_long$Percentage), "%"
)

# 7. Create the figure
p <- ggplot(
  plot_long,
  aes(x = X, y = Percentage, fill = interaction(Class, Category))
) +
  
  geom_col(
    width = 0.42,
    color = "white",
    linewidth = 0.35
  ) +
  
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "EAA.RQAA > 1" = "#B33B4A",
      "EAA.RQAA \u2264 1" = "#F0C7CC",
      "non-EAA.RQAA > 1" = "#1769AA",
      "non-EAA.RQAA \u2264 1" = "#BFD8F0"
    ),
    breaks = c(
      "EAA.RQAA > 1",
      "EAA.RQAA \u2264 1",
      "non-EAA.RQAA > 1",
      "non-EAA.RQAA \u2264 1"
    ),
    labels = c(
      "EAA: RQAA > 1",
      "EAA: RQAA \u2264 1",
      "non-EAA: RQAA > 1",
      "non-EAA: RQAA \u2264 1"
    ),
    name = "Relative quantity"
  ) +
  
  scale_x_continuous(
    breaks = group_centers,
    labels = group_order,
    limits = c(0.4, 8.1),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    title = "Relative amino-acid quantities across taxonomic groups",
    x = NULL,
    y = "Measurements (%)"
  ) +
  
  # Add EAA/non-EAA labels beneath each bar
  annotate(
    "text",
    x = c(
      group_centers[1] - 0.23, group_centers[1] + 0.23,
      group_centers[2] - 0.23, group_centers[2] + 0.23,
      group_centers[3] - 0.23, group_centers[3] + 0.23,
      group_centers[4] - 0.23, group_centers[4] + 0.23,
      group_centers[5] - 0.23, group_centers[5] + 0.23
    ),
    y = -7,
    label = rep(c("EAA", "non-EAA"), 5),
    size = 3.2
  ) +
  
  coord_cartesian(
    ylim = c(0, 100),
    clip = "off"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5,
     margin = margin(b = 20)
    ),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 10
    ),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(0.55, "cm"),
    plot.margin = margin(10, 10, 30, 10)
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

# 8. Display the graph
print(p)

# 9. Save high-resolution PNG and PDF
ggsave(
  "RQAA_taxonomic_groups.png",
  plot = p,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "RQAA_taxonomic_groups.pdf",
  plot = p,
  width = 6,
  height = 6.5,
  units = "in",
  device = "pdf",
  bg = "white"
)



