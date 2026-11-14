
# FIGURA 1 — RQAA por grupo taxonomico
# Puntos individuales + mediana + intervalo intercuartilico

library(tidyverse)
library(readxl)
library(janitor)

# ---- CONFIGURACION ----
archivo_datos <- "dataFigure1.csv"
hoja_excel <- 1
carpeta_salida <- "Your Output Path Here/"

# ---- IMPORTAR DATOS ----
if (grepl("\\.xlsx?$", archivo_datos, ignore.case = TRUE)) {
  datos_originales <- read_excel(archivo_datos, sheet = hoja_excel)
} else if (grepl("\\.csv$", archivo_datos, ignore.case = TRUE)) {
  datos_originales <- read_csv(
    archivo_datos,
    na = c("", "NA"),
    show_col_types = FALSE
  )
} else {
  stop("Usa un archivo Excel o CSV.")
}

datos <- datos_originales |> clean_names()

required <- c("group", "amino_acid", "rqaa", "class")
missing_cols <- setdiff(required, names(datos))

if (length(missing_cols) > 0) {
  stop(paste(
    "Faltan columnas:", paste(missing_cols, collapse = ", "),
    ". Se requieren group, amino_acid, rqaa y class."
  ))
}

# ---- LIMPIEZA ----
datos <- datos |>
  mutate(
    group = str_squish(as.character(group)),
    amino_acid = str_to_title(str_squish(as.character(amino_acid))),
    class = str_to_lower(str_squish(as.character(class))),
    rqaa = suppressWarnings(as.numeric(rqaa)),
    class = case_when(
      class %in% c(
        "eaa", "essential", "essential amino acid",
        "essential amino acids"
      ) ~ "EAA",
      class %in% c(
        "non-eaa", "non eaa", "non-essential",
        "non essential", "non-essential amino acid",
        "non-essential amino acids"
      ) ~ "Non-EAA",
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(group), group != "",
    !is.na(amino_acid), amino_acid != "",
    !is.na(rqaa), is.finite(rqaa),
    !is.na(class)
  )

if (nrow(datos) == 0) {
  stop("No quedan filas validas despues de limpiar los datos.")
}

# Los ceros representan datos ausentes en este conjunto.
datos <- datos |>
  mutate(rqaa = if_else(rqaa == 0, NA_real_, rqaa)) |>
  filter(!is.na(rqaa), rqaa > 0)

# ---- ORDEN DE GRUPOS Y AMINOACIDOS ----
group_order <- c(
  "Plants","Bacteria", "Fungi", "Animals"
)

# Orden de clase: no esenciales primero y esenciales despues.
orden_aa <- c(
  "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu",
  "Gly", "Pro", "Ser", "Tyr",
  "His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val"
)

orden_por_clase <- datos |>
  distinct(amino_acid, class) |>
  mutate(
    aa_texto = as.character(amino_acid),
    clase_num = match(class, c("Non-EAA", "EAA")),
    orden_num = match(aa_texto, orden_aa)
  ) |>
  arrange(clase_num, orden_num, aa_texto) |>
  pull(aa_texto) |>
  unique()

datos_plot <- datos |>
  mutate(
    amino_acid = factor(amino_acid, levels = orden_por_clase),
    group = factor(group, levels = group_order),
    class = factor(class, levels = c("Non-EAA", "EAA"))
  )

# ---- RESUMEN POR GRUPO Y AMINOACIDO ----
resumen_aa <- datos_plot |>
  group_by(group, amino_acid, class) |>
  summarise(
    n = sum(is.finite(rqaa) & rqaa > 0),
    mediana = median(rqaa[rqaa > 0 & is.finite(rqaa)]),
    q25 = unname(quantile(rqaa[rqaa > 0 & is.finite(rqaa)], 0.25)),
    q75 = unname(quantile(rqaa[rqaa > 0 & is.finite(rqaa)], 0.75)),
    .groups = "drop"
  ) |>
  filter(n > 0)

n_no_eaa <- datos_plot |>
  distinct(amino_acid, class) |>
  filter(class == "Non-EAA") |>
  nrow()

n_eaa <- datos_plot |>
  distinct(amino_acid, class) |>
  filter(class == "EAA") |>
  nrow()

# ---- RANGO DEL EJE Y ----
rango_y <- range(datos_plot$rqaa, na.rm = TRUE)
y_min <- min(0.01, 10^floor(log10(rango_y[1])))
y_max <- max(10, 10^ceiling(log10(rango_y[2])))

# ---- CONSTRUIR FIGURA ----
figura1 <- ggplot(
  datos_plot,
  aes(x = amino_acid, y = rqaa, color = class)
) +
  # Fondo de clases
  annotate(
    "rect",
    xmin = 0.5, xmax = n_no_eaa + 0.5,
    ymin = y_min, ymax = y_max,
    fill = "#E8F1FF", alpha = 0.9
  ) +
  annotate(
    "rect",
    xmin = n_no_eaa + 0.5, xmax = n_no_eaa + n_eaa + 0.5,
    ymin = y_min, ymax = y_max,
    fill = "#FCE8E8", alpha = 0.9
  ) +
  # Referencia RQAA = 1
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1,
    label = "RQAA = 1",
    hjust = 1.1,
    vjust = -0.5,
    size = 4.5,
    fontface = "bold",
    color = "grey25"
  ) +
  # Datos individuales
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 123),
    alpha = 0.40,
    size = 1.1
  ) +
  
  # Intervalo intercuartilico: linea fina con extremos
  geom_errorbar(
    data = resumen_aa,
    aes(
      x = amino_acid,
      ymin = q25,
      ymax = q75,
      color = class
    ),
    inherit.aes = FALSE,
    width = 0.16,
    linewidth = 0.55,
    alpha = 1
  ) +
  # Mediana: punto pequeño e intenso del color de cada clase
  geom_point(
    data = resumen_aa,
    aes(
      x = amino_acid,
      y = mediana,
      color = class
    ),
    inherit.aes = FALSE,
    size = 1.1,
    alpha = 1
  ) +
  facet_wrap(
    ~ group,
    ncol = 2,
    drop = TRUE,
    axes = "all",
    axis.labels = "all"
  ) +
  
  scale_color_manual(
    values = c(
      "Non-EAA" = "#0066A8",
      "EAA" = "#C62828"
    ),
    breaks = c("Non-EAA", "EAA"),
    labels = c(
      "Non-essential amino acids",
      "Essential amino acids"
    ),
    name = NULL
  ) +
  
  scale_y_log10(
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = scales::label_number(),
    limits = c(y_min, y_max)
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = "RQAA (log scale)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 13, color = "grey20"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 13,
      color = "grey15"
    ),
    legend.position = "none",
    panel.spacing = grid::unit(1.5, "lines"),
    plot.margin = margin(8, 10, 8, 8)
  )

print(figura1)

# ---- EXPORTAR ----

dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(carpeta_salida, "Figure1_RQAA.tiff"),
  plot = figura1,
  device = "tiff",
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write_csv(
  datos_plot |>
    mutate(
      group = as.character(group),
      amino_acid = as.character(amino_acid),
      class = as.character(class)
    ),
  file.path(carpeta_salida, "Figure1_data_cleaned.csv")
)

cat("\nFigura generada. Revisa la carpeta de salida y los avisos.\n")


