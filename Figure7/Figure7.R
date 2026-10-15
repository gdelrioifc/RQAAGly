library(scatterplot3d)

setwd("/Volumes/Disco2024/Papers/EnPreparacion/DProtein/2026/EAAInProteinsAndMetabolism/DataFiguresPaper/Figure7/OurWorldInData")

data1 <- read.delim("Figure7A.tsv", header = TRUE, sep = "\t")
data2 <- read.delim("Figure7B.tsv", header = TRUE, sep = "\t")

# Determine axis limits dynamically so both datasets fit comfortably within the plot boundaries
x_range <- range(c(data1$GHE, data2$GHE))
y_range <- range(c(data1$LandUse, data2$LandUse))
z_range <- range(c(data1$WaterUse, data2$WaterUse))

# Create the initial plot using the first dataset
s3d <- scatterplot3d(data1$GHE, data1$LandUse, data1$WaterUse,
                     main = "",
                     xlab = "GHE",
                     ylab = "Land Use",
                     zlab = "Water Use",
                     xlim = x_range,
                     ylim = y_range,
                     zlim = z_range,
                     color = "dodgerblue", # Color for File 1
                     pch = 19,
                     type = "p",
                     z.ticklabs = rep("", length(z_ticks))) # Hides Z values completely

# Use s3d$points3d to overlay the second dataset onto the same plot in a different color
s3d$points3d(data2$GHE, data2$LandUse, data2$WaterUse, 
             color = "firebrick",  # Color for File 2
             pch = 17)             # Different point shape (triangle) for visual distinction


# Add a legend to clarify which color belongs to which file
# legend("topright", 
#       legend = c("File 1 (Group A)", "File 2 (Group B)"),
#       col = c("dodgerblue", "firebrick"), 
#       pch = c(19, 17), 
#       bg = "white")