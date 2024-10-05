#set working directory
setwd("C:/Users/PMYLS/Desktop/Assignment")
#loading files
library(readxl)
RAW_DATA_ <- read_excel("RAW DATA .xlsx")
#Genetic markers start from column 6 to the last column
genetic_data <- RAW_DATA_[, 6:ncol(RAW_DATA_)]

# Function to calculate Shannon's information index (I)
shannon_info <- function(x) {
  -sum(table(x)/length(x) * log(table(x)/length(x)), na.rm = TRUE)
}

# Calculate key descriptive statistics
shannon_index <- apply(genetic_data, 1, shannon_info)
obs_heterozygosity <- apply(genetic_data, 1, function(x) sum(x != "A_A") / length(x))
exp_heterozygosity <- apply(genetic_data, 1, function(x) 2 * sum(x == "A_A") * sum(x == "A_T") / length(x)^2)
minor_allele_freq <- apply(genetic_data, 1, function(x) sum(x == "A_T") / length(x))
fixation_index <- apply(genetic_data, 1, function(x) (sum(x == "A_T") / length(x)) - (sum(x == "A_A") / length(x)))
inbreeding <- apply(genetic_data, 1, function(x) 1 - obs_heterozygosity)

# Create a data frame with the calculated statistics
statistics_df <- data.frame(
  Accessions = RAW_DATA_$Accessions,
  Province = RAW_DATA_$Province,
  Shannon_Info = shannon_index,
  Obs_Heterozygosity = obs_heterozygosity,
  Exp_Heterozygosity = exp_heterozygosity,
  Minor_Allele_Freq = minor_allele_freq,
  Fixation_Index = fixation_index,
  Inbreeding = inbreeding
)
write.csv(statistics_df, "statisics_df.csv")

# Install and load required packages
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan")
}
if (!requireNamespace("geometry", quietly = TRUE)) {
  install.packages("geometry")
}
library(vegan)
library(geometry)

# Extract the relevant columns for analysis
data_for_pcoa <- statistics_df[, c("Accessions", "Shannon_Info", "Obs_Heterozygosity", "Exp_Heterozygosity", "Minor_Allele_Freq", "Fixation_Index")]

# Set the 'Accessions' column as row names
row.names(data_for_pcoa) <- data_for_pcoa$Accessions
data_for_pcoa <- data_for_pcoa[, -1]  # Remove the 'Accessions' column

# Perform PCoA
dist_matrix <- vegdist(data_for_pcoa, method = "euclidean")
pcoa_result <- cmdscale(dist_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

# Cut the tree into k clusters (adjust k as needed)
k_clusters <- cutree(hc, k = 3)

# Plot dendrogram
plot(hc, main = "Dendrogram of Hierarchical Clustering")

# Plot PCoA with colored clusters
colors <- c("red", "green", "blue")  # Adjust as needed
plot(pcoa_result[, 1], pcoa_result[, 2], main = "Principal Coordinate Analysis with Clusters", 
     xlab = "PCo1", ylab = "PCo2", pch = 19, col = colors[k_clusters])

# Add labels for each point
text(pcoa_result[, 1], pcoa_result[, 2], labels = row.names(data_for_pcoa), pos = 3, col = "black")

# Encircle clusters and label them
for (i in unique(k_clusters)) {
  cluster_points <- pcoa_result[k_clusters == i, c(1, 2)]
  hull <- convhulln(cluster_points)
  polygon(hull, col = NA, border = colors[i], lwd = 2)
  center <- colMeans(cluster_points)
  text(center[1], center[2], labels = paste("Cluster", i), col = colors[i], pos = 1, font = 2)
}

# Print the eigenvalues
summary(pcoa_result)


# Install and load required packages
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(pheatmap)
library(RColorBrewer)

#Extract relevant columns
heatmap_data <- statistics_df[, c("Shannon_Info", "Obs_Heterozygosity", "Exp_Heterozygosity", "Minor_Allele_Freq", "Fixation_Index")]

# Create a heatmap with pheatmap
pheatmap(heatmap_data, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of Variable Distributions",
         annotation_names_col = TRUE)


if (!requireNamespace("GGally", quietly = TRUE)) {
  install.packages("GGally")
}

library(GGally)

# Extract relevant columns
scatter_data <- statistics_df[, c("Shannon_Info", "Obs_Heterozygosity", "Exp_Heterozygosity", "Minor_Allele_Freq", "Fixation_Index")]

# Add colors based on a categorical variable (e.g., cluster information)
scatter_data$Cluster <- as.factor(k_clusters)

# Create a scatterplot matrix
ggpairs(
  scatter_data,
  mapping = aes(color = Cluster),
  title = "Scatterplot Matrix of Variables",
  lower = list(continuous = wrap("points", alpha = 0.5)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
)

# Subset the data for each province
bocas_del_toro_data <- subset(statistics_df, Province == "Bocas del Toro")
herrera_data <- subset(statistics_df, Province == "Herrera")
panama_data <- subset(statistics_df, Province == "Panama")
cocle_data <- subset(statistics_df, Province == "Coclé")
colon_data <- subset(statistics_df, Province == "Colón")
chiriqui_data <- subset(statistics_df, Province == "Chiriquí")
darien_data <- subset(statistics_df, Province == "Darién")
veraguas_data <- subset(statistics_df, Province == "Veraguas")

# Summary statistics
summary(bocas_del_toro_data)
summary(herrera_data)
summary(panama_data)
summary(cocle_data)
summary(colon_data)
summary(chiriqui_data)
summary(darien_data)
summary(veraguas_data)

# Boxplot for Shannon_Info across provinces
ggplot(statistics_df, aes(x = Province, y = Shannon_Info)) +
  geom_boxplot() +
  labs(title = "Boxplot of Shannon_Info by Province")

# Scatter plot for Obs_Heterozygosity and Exp_Heterozygosity
ggplot(statistics_df, aes(x = Obs_Heterozygosity, y = Exp_Heterozygosity, color = Province)) +
  geom_point() +
  labs(title = "Scatter Plot of Obs_Heterozygosity vs. Exp_Heterozygosity by Province")

# Pairwise scatter plots for select variables
selected_vars <- c("Shannon_Info", "Obs_Heterozygosity", "Exp_Heterozygosity", "Minor_Allele_Freq")
pairs(statistics_df[selected_vars], main = "Pairwise Scatter Plots of Select Variables")

# Load required libraries
library(circlize)

# Extract unique provinces and create a table
province_table <- table(statistics_df$Province)

# Convert table to a data frame
province_counts <- data.frame(province = names(province_table), count = as.numeric(province_table))

# Sort by province count
province_counts <- province_counts[order(province_counts$count, decreasing = TRUE), ]

# Create a matrix with hierarchical structure
province_matrix <- matrix(0, nrow = length(unique_provinces), ncol = length(unique_provinces))
rownames(province_matrix) <- colnames(province_matrix) <- unique_provinces

# Fill the matrix with province counts
for (i in 1:nrow(province_matrix)) {
  province_matrix[i, ] <- province_counts$count
}

# Create a circle cluster tree
chordDiagram(province_matrix, transparency = 0.5)
# Density plot for Fixation Index by province
library(ggplot2)
ggplot(statistics_df, aes(x = Fixation_Index, fill = Province)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Fixation Index by Province", x = "Fixation Index")
library(tidyr)
library(dplyr)


# Select inbreeding columns along with the "Province" column
inbreeding_data <- select(statistics_df, Province, starts_with("Inbreeding"))

# Gather the data
inbreeding_data_long <- gather(inbreeding_data, key = "Inbreeding_Coefficient", value = "Value")

# Add a third column with the province for each inbreeding coefficient
inbreeding_data_long$Province_Name <- statistics_df$Province[rep(seq_len(nrow(statistics_df)), each = ncol(inbreeding_data))]

# View the resulting data frame
head(inbreeding_data_long)

# Remove the first 153 rows
inbreeding_data_long <- inbreeding_data_long[-(1:153), ]

# View the modified data frame
head(inbreeding_data_long)
# Stacked bar plot for accessions count by province and inbreeding coefficient
ggplot(inbreeding_data_long, aes(x = Inbreeding_Coefficient, fill = Province_Name)) +
  geom_bar() +
  labs(title = "Stacked Bar Plot of Accessions Count by Inbreeding Coefficient and Province", x = "Inbreeding Coefficient", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Heatmap for inbreeding coefficients
ggplot(inbreeding_data_long, aes(x = Inbreeding_Coefficient, y = Province_Name )) +
  geom_tile(aes(fill = Province_Name)) +
  labs(title = "Heatmap of Inbreeding Coefficients by Province", x = "Province", y = "Inbreeding Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
