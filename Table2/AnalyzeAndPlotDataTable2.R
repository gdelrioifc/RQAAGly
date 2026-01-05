library(corrplot)
setwd([The path where to find the RDA_RQAA CSV file])
# Please use the sheets labeled as RDA_RQAAGly/Val/Leu in file SupplementalData_Table2.xlsx
data <- read.csv("RDA_RQAAGly.csv")
correlation_matrix <- cor(data, use = "pairwise.complete.obs")
# Save the correlation matrix as a CSV file
write.csv(correlation_matrix, "RDA_RQAAGly_correlation_matrix.csv", row.names = TRUE) 
# display correlation matrix as a heatmap
heatmap(correlation_matrix)
