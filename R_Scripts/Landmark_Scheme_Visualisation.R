# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Landmark Scheme Visualisation 

################################################################################

# Load the necessary libraries 
library(Morpho)
library(abind)
library(rgl)
library(Rvcg)
library(dplyr)
library(geomorph)

#########################################################################################
# Load Comparison Meshes

# Load Comparison Mesh 1 
mesh1 <- read.ply("E:/James Mulqueeney/Paper 2- Mammal Shape/Arctictis_binturong_Decimated.ply")

#########################################################################################

# Manual Landmarking Heatmaps

# Read in the landmark data 
shape_data <- read.csv("E:/James Mulqueeney/Paper 2- Mammal Shape/Write Up/BMC Ecology and Evolution/Data & Code/Data/Data/Data_A18-Mirrored_322.csv", header = TRUE)

# Extract x, y, z coordinates in 3D
x_values <- as.matrix(shape_data[, grep("\\.X", colnames(shape_data))])
y_values <- as.matrix(shape_data[, grep("\\.Y", colnames(shape_data))])
z_values <- as.matrix(shape_data[, grep("\\.Z", colnames(shape_data))])

num_landmarks <- ncol(x_values)
num_specimens <- nrow(shape_data)

# Initialize result matrix
result_matrix <- array(0, dim = c(num_landmarks, 3, num_specimens))

# Populate result matrix with x, y, z values
for (i in 1:num_specimens) {
  result_matrix[, 1, i] <- x_values[i, ]
  result_matrix[, 2, i] <- y_values[i, ]
  result_matrix[, 3, i] <- z_values[i, ]
}

# Concatenate x, y, z values
shape_concatenated_array <- abind(result_matrix, along = 1)

# Plot the atlas data
spheres3d(shape_concatenated_array[,,22], radius = 1) 

#########################################################################################

# Visualisation 

#open a 3d window with two viewing frames
open3d()

# Plot the mesh 
shade3d(mesh1,col="#DECAB0")

num_rows <- nrow(shape_concatenated_array[,,22])  # Total number of points

# Define color vector: First 60 in red, 755-820 in red, rest in gold
colors <- rep("gold", num_rows)  # Default all to gold
colors[1:66] <- "red"            # First 60 points in red
colors[755:(754+51)] <- "red"    # 66 points after 754 in red

# Plot the points with the defined colors
spheres3d(shape_concatenated_array[,,22], radius = 1, col = colors)

#########################################################################################
