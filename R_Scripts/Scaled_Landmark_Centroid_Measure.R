# Load necessary libraries
library(plotly)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)
library(geomorph)
library(ape)
library(Morpho)

#########################################################################################

# Manual Landmarking Results 

#########################################################################################

# Read in the CSV file
data <- read.csv("C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/Data_A24-Combined_Data_322.csv", stringsAsFactors = FALSE)

# Re-organise the shape matrix to fit geomorph specifications 

# Extract the shape coordinates (assuming columns ending with '.X', '.Y', or '.Z')
coordinate_columns <- grep("\\.X|\\.Y|\\.Z", names(data), value = TRUE)

# Convert the coordinate columns to a matrix
shape_data_matrix <- as.matrix(data[, coordinate_columns])

# Determine the number of specimens and landmarks
num_specimens <- nrow(shape_data_matrix)  # Number of rows (specimens)
num_dimensions <- 3  # Typically 3 for X, Y, Z
num_landmarks <- length(coordinate_columns) / num_dimensions  # Calculate number of landmarks

# Reshape the matrix into a 3D array: Specimens × Landmarks × Dimensions
shape_data_3d <- arrayspecs(shape_data_matrix, p = num_landmarks, k = num_dimensions)

# Assign the tip labels to the matrix 
dimnames(shape_data_3d) <- list(NULL,  NULL, as.character(data$Tip_Label)) 

# Construct the geomorph data frame with appropriate checks
gdf <- geomorph.data.frame(
  coords = shape_data_3d,
  specimen_id = as.factor(data$Tip_Label),  # Ensure factor type for grouping variables
  diet = as.factor(data$Diet), 
  loc = as.factor(data$Locomotion), 
  dev = as.factor(data$Development), 
  clade = as.factor(data$Order), 
  superorder = as.factor(data$Superorder)
)

# Verify the structure of the geomorph data frame
print(gdf)
str(gdf)

#########################################################################################

# Calculate centroids 

# Initialize a vector to store centroid sizes
centroid_sizes <- numeric(dim(gdf$coords)[3])

# Loop through each specimen and calculate centroid size
for (i in 1:dim(gdf$coords)[3]) {
  centroid_sizes[i] <- cSize(gdf$coords[,,i])
}

# Create a data frame to save the results
results <- data.frame(Specimen = 1:length(centroid_sizes), CentroidSize = centroid_sizes)

# Save the results to \a CSV file
write.csv(results, "E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Deterministic-Atlas-Analysis-main (2)/Deterministic-Atlas-Analysis-main/Data/New Data/Centroid Data/Centroid_Sizes2.csv", row.names = FALSE)
