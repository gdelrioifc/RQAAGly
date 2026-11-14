
# ============================================================
# FIGURE 3. PROTEOMIC RQAA DISTRIBUTION ACROSS 23 SPECIES
# ============================================================

library(tidyverse)

# ============================================================
# 1. FILE PATHS
# ============================================================

input_file <- "dataFigure3.csv"

output_dir <- "Your Ouput Path Here"

# ============================================================
# 2. IMPORT DATA
# ============================================================

data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  fileEncoding = "Windows-1252"
)

# ============================================================
# 3. DEFINE AMINO ACID CLASSES
# ============================================================

essential_AA <- c(
  "Phe", "His", "Ile", "Leu", "Lys",
  "Met", "Thr", "Trp", "Val"
)

nonessential_AA <- c(
  "Ala", "Asp", "Glu", "Asn", "Pro", "Gln",
  "Arg", "Ser", "Gly", "Cys", "Tyr"
)

AA_order <- c(nonessential_AA, essential_AA)

# ============================================================
# 4. CONVERT TO LONG FORMAT
# ============================================================

data_long <- data %>%
  pivot_longer(
    cols = -c(Index, AARatio),
    names_to = "Species",
    values_to = "RQAA"
  ) %>%
  rename(AminoAcid = AARatio) %>%
  mutate(
    RQAA = as.numeric(RQAA),
    
    AminoAcid = factor(
      AminoAcid,
      levels = AA_order
    ),
    
    Class = case_when(
      AminoAcid %in% essential_AA ~ "Essential amino acids",
      AminoAcid %in% nonessential_AA ~ "Non-essential amino acids"
    ),
    
    Class = factor(
      Class,
      levels = c(
        "Non-essential amino acids",
        "Essential amino acids"
      )
    )
  ) %>%
  filter(!is.na(RQAA))

# Check data
cat("Number of species:", n_distinct(data_long$Species), "\n")
cat("Number of measurements:", nrow(data_long), "\n")

# ============================================================
# 5. CALCULATE MEDIAN AND IQR
# ============================================================

data_summary <- data_long %>%
  group_by(AminoAcid, Class) %>%
  summarise(
    Median = median(RQAA, na.rm = TRUE),
    Q1 = quantile(RQAA, 0.25, na.rm = TRUE),
    Q3 = quantile(RQAA, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 6. BACKGROUND COLORS
# ============================================================

background <- tibble(
  xmin = c(0.5, 11.5),
  xmax = c(11.5, 20.5),
  fill = c("#DCEBFA", "#FBE0E0")
)

# ============================================================
# 7. GRAPH
# ============================================================

p <- ggplot() +
  
  # Background colors
  geom_rect(
    data = background,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = background$fill,
    inherit.aes = FALSE
  ) +
  
  # Reference line: Glycine = 1
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.55
  ) +
  
  # Individual proteomes
  geom_jitter(
    data = data_long,
    aes(
      x = AminoAcid,
      y = RQAA,
      color = Class,
      shape = "Individual proteomes"
    ),
    width = 0.12,
    height = 0,
    size = 1.35,
    alpha = 0.55
  ) +
  
  # Interquartile range (Q1-Q3)
  geom_errorbar(
    data = data_summary,
    aes(
      x = AminoAcid,
      ymin = Q1,
      ymax = Q3,
      color = Class
    ),
    width = 0.08,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  
  # Median 
  geom_point(
    data = data_summary,
    aes(
      x = AminoAcid,
      y = Median,
      color = Class
    ),
    shape = 16,
    size = 2.2
  ) +
  
  
  # Color legend: amino acid class
  scale_color_manual(
    name = NULL,
    values = c(
      "Non-essential amino acids" = "#174A9C",
      "Essential amino acids" = "#C71919"
    ),
    breaks = c(
      "Non-essential amino acids",
      "Essential amino acids"
    ),
    guide = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        size = 3,
        alpha = 1
      )
    )
  ) +
  
  
  # X axis
  scale_x_discrete(
    limits = AA_order,
    drop = FALSE
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = "RQAA (proteome)"
  ) +
  
  guides(
    color = "none",
    shape = "none"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 15,
      face = "bold",
      margin = margin(r = 8)
    ),
    
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.45
    ),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.45
    ),
    
    legend.position = "none",
    
    plot.margin = margin(10, 12, 8, 12)
  )

# 9. DISPLAY GRAPH

print(p)

# ============================================================
# 10. EXPORT TIFF
# ============================================================

ggsave(
  filename = file.path(
    output_dir,
    "Figure3_Proteome_RQAA.tiff"
  ),
  plot = p,
  width = 10,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 11. EXPORT PDF
# ============================================================

ggsave(
  filename = file.path(
    output_dir,
    "Figure3_Proteome_RQAA.pdf"
  ),
  plot = p,
  width = 10,
  height = 5.5,
  bg = "white"
)
