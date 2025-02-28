# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Manual landmarking heatmaps 

################################################################################

# Load required libraries 
library(Morpho)
library(abind)
library(rgl)
library(Rvcg)
library(dplyr)
library(geomorph)

################################################################################

# Manual Landmarking Heatmaps

################################################################################

# Read in the landmark data 
shape_data <- read.csv("C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/Data_A21-Shape_Data_322.csv", header = TRUE)

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

# Calculate the average for each column across all matrix slices
shape_average_matrix_slice <- apply(shape_concatenated_array, MARGIN = c(1, 2), FUN = mean)

# Plot the average matrix slice
spheres3d(shape_average_matrix_slice, radius = 0.001, col = "Red")

################################################################################

# Generate and plot mirrored coordinates for the mean shape 

# Mirror the x-values
shape_average_matrix_slice_mirrored <- shape_average_matrix_slice
shape_average_matrix_slice_mirrored[, 1] <- -shape_average_matrix_slice[, 1] 

# Rotate around YZ plane (180 degrees rotation around X-axis)
shape_average_matrix_slice_mirrored[, 2] <- -shape_average_matrix_slice_mirrored[, 2]
shape_average_matrix_slice_mirrored[, 1] <- -shape_average_matrix_slice_mirrored[, 1]

# Plot the mirrored shape
spheres3d(shape_average_matrix_slice_mirrored, radius = 0.001, col = "Blue")

# Combine the mirrored and original matrices
combined_matrix0 <- rbind(shape_average_matrix_slice, shape_average_matrix_slice_mirrored)
spheres3d(combined_matrix0, radius = 0.001)

################################################################################

# Generate mirrored coordinates for all shapes

# Create a copy of the original array to store the transformed data
shape_concatenated_array2_transformed <- shape_concatenated_array

# Iterate through each slice of the array
for (i in 1:dim(shape_concatenated_array)[3]) {
  # Extract the current slice
  current_slice <- shape_concatenated_array[,,i]
  
  # Apply transformations
  current_slice[, 1] <- -current_slice[, 1]  # Negate X-values
  current_slice[, 2] <- -current_slice[, 2]  # Negate Y-values
  current_slice[, 1] <- -current_slice[, 1]  # Rotate around YZ plane (180 degrees rotation around X-axis)
  
  # Save the transformed slice back into the transformed array
  shape_concatenated_array2_transformed[,,i] <- current_slice
}

################################################################################

# Define input directory
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Outputs/"

################################################################################

# Apply to Artiodactyla Specimen - Bos domestic

################################################################################

# Load Comparison Mesh
mesh1 <- read.ply(file.path(input_dir, "Bos_domestic.ply"))

# Combine matrices of the selected specimen
combined_matrix1 <- rbind(shape_concatenated_array[,,32], shape_concatenated_array2_transformed[,,32])
spheres3d(combined_matrix1, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh1_warped <-plotRefToTarget(combined_matrix1, combined_matrix0, mesh = mesh1, method = "surface" )

# Compare the two meshes
shade3d(mesh1_warped, col = 2, specular = 1)
shade3d(mesh1, col = 3, specular = 1)

# Visualize within the heatmaps
MD1 <- meshDist(mesh1, mesh1_warped, from = -150, to = 20)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh1_warped <- tps3d(mesh1, combined_matrix1, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh1_warped, col = 2, specular = 1)
shade3d(mesh1, col = 3, specular = 1)

# Visualize within the heatmaps
MD1 <- meshDist(mesh1, mesh1_warped, from = -150, to = 20)

#########################################################################################

# Apply to Carnivora - Arctictis binturong

#########################################################################################

# Load Comparison Mesh 
mesh2 <- read.ply(file.path(input_dir, "Arctictis_binturong.ply"))

# Combine matrices of the selected specimen
combined_matrix2 <- rbind(shape_concatenated_array[,,22], shape_concatenated_array2_transformed[,,22])
spheres3d(combined_matrix2, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh2_warped <-plotRefToTarget(combined_matrix2, combined_matrix0, mesh = mesh2, method = "surface" )

# Compare the two meshes
shade3d(mesh2_warped, col = 2, specular = 1)
shade3d(mesh2, col = 3, specular = 1)

# Visualize within the heatmaps
MD2 <- meshDist(mesh2, mesh2_warped, from = -5, to = 5)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh2_warped <- tps3d(mesh2, combined_matrix2, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh2_warped, col = 2, specular = 1)
shade3d(mesh2, col = 3, specular = 1)

# Visualize within the heatmaps
MD2 <- meshDist(mesh2, mesh2_warped, from = -5, to = 5)

#########################################################################################

# Apply to Cetacea - Schizodelphis morckhoviensis

#########################################################################################

# Load Comparison Mesh  
mesh3 <- read.ply(file.path(input_dir, "Schizodelphis_morckhoviensis.ply"))

# Combine matrices of the selected specimen
combined_matrix3 <- rbind(shape_concatenated_array[,,272], shape_concatenated_array2_transformed[,,272])
spheres3d(combined_matrix3, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh3_warped <-plotRefToTarget(combined_matrix3, combined_matrix0, mesh = mesh3, method = "surface" )

# Compare the two meshes
shade3d(mesh3_warped, col = 2, specular = 1)
shade3d(mesh3, col = 3, specular = 1)

# Visualize within the heatmaps
MD2 <- meshDist(mesh3, mesh3_warped, from=-80, to = 20)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh3_warped <- tps3d(mesh3, combined_matrix3, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh3_warped, col = 2, specular = 1)
shade3d(mesh3, col = 3, specular = 1)

# Visualize within the heatmaps
MD3 <- meshDist(mesh3, mesh3_warped, from=-80, to = 20)

#########################################################################################

# Apply to Chiroptera - Centurio senex

#########################################################################################

# Load Comparison Mesh  
mesh4 <- read.ply(file.path(input_dir, "Centurio_senex.ply"))

# Combine matrices of the selected specimen
combined_matrix4 <- rbind(shape_concatenated_array[,,54], shape_concatenated_array2_transformed[,,54])
spheres3d(combined_matrix4, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh4_warped <-plotRefToTarget(combined_matrix4, combined_matrix0, mesh = mesh4, method = "surface" )

# Compare the two meshes
shade3d(mesh4_warped, col = 2, specular = 1)
shade3d(mesh4, col = 3, specular = 1)

# Visualize within the heatmaps
MD4 <- meshDist(mesh4, mesh4_warped, from=-40, to = 15)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh4_warped <- tps3d(mesh4, combined_matrix4, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh4_warped, col = 2, specular = 1)
shade3d(mesh4, col = 3, specular = 1)

# Visualize within the heatmaps
MD4 <- meshDist(mesh4, mesh4_warped, from=-40, to = 15)

#########################################################################################

# Apply to Perrisodactyla - Equus caballus

#########################################################################################

# Load Comparison Mesh  
mesh5 <- read.ply(file.path(input_dir, "Equus_caballus.ply"))

# Combine matrices of the selected specimen
combined_matrix5 <- rbind(shape_concatenated_array[,,108], shape_concatenated_array2_transformed[,,108])
spheres3d(combined_matrix5, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh5_warped <-plotRefToTarget(combined_matrix5, combined_matrix0, mesh = mesh5, method = "surface" )

# Compare the two meshes
shade3d(mesh5_warped, col = 2, specular = 1)
shade3d(mesh5, col = 3, specular = 1)

# Visualize within the heatmaps
MD5 <- meshDist(mesh5, mesh5_warped, from=-40, to = 15)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh5_warped <- tps3d(mesh5, combined_matrix5, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh5_warped, col = 2, specular = 1)
shade3d(mesh5, col = 3, specular = 1)

# Visualize within the heatmaps
MD5 <- meshDist(mesh5, mesh5_warped, from=-40, to = 15)

#########################################################################################

# Apply to Pilosa - Choloepus hoffmanni

#########################################################################################

# Load Comparison Mesh  
mesh6 <- read.ply(file.path(input_dir, "Choloepus_hoffmanni.ply"))

# Combine matrices of the selected specimen
combined_matrix6 <- rbind(shape_concatenated_array[,,61], shape_concatenated_array2_transformed[,,61])
spheres3d(combined_matrix6, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh6_warped <-plotRefToTarget(combined_matrix6, combined_matrix0, mesh = mesh6, method = "surface" )

# Compare the two meshes
shade3d(mesh6_warped, col = 2, specular = 1)
shade3d(mesh6, col = 3, specular = 1)

# Visualize within the heatmaps
MD6 <- meshDist(mesh6, mesh6_warped, from=-80, to = 20)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh6_warped <- tps3d(mesh6, combined_matrix6, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh6_warped, col = 2, specular = 1)
shade3d(mesh6, col = 3, specular = 1)

# Visualize within the heatmaps
MD6 <- meshDist(mesh6, mesh6_warped, from=-80, to = 20)

#########################################################################################

# Apply to Primates - Cacajao_calvus

#########################################################################################

# Load Comparison Mesh  
mesh7 <- read.ply(file.path(input_dir, "Cacajao_calvus.ply"))

# Combine matrices of the selected specimen
combined_matrix7 <- rbind(shape_concatenated_array[,,36], shape_concatenated_array2_transformed[,,36])
spheres3d(combined_matrix7, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh7_warped <-plotRefToTarget(combined_matrix7, combined_matrix0, mesh = mesh7, method = "surface" )

# Compare the two meshes
shade3d(mesh7_warped, col = 2, specular = 1)
shade3d(mesh7, col = 3, specular = 1)

# Visualize within the heatmaps
MD7 <- meshDist(mesh7, mesh7_warped, from=-40, to = 15)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh7_warped <- tps3d(mesh7, combined_matrix7, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh7_warped, col = 2, specular = 1)
shade3d(mesh7, col = 3, specular = 1)

# Visualize within the heatmaps
MD7 <- meshDist(mesh7, mesh7_warped, from=-40, to = 15)

#########################################################################################

# Apply to Rodentia - Cavia aperea

#########################################################################################

# Load Comparison Mesh  
mesh8 <- read.ply(file.path(input_dir, "Cavia_aperea.ply"))

# Combine matrices of the selected specimen
combined_matrix8 <- rbind(shape_concatenated_array[,,50], shape_concatenated_array2_transformed[,,51])
spheres3d(combined_matrix8, radius = 0.001, col = "Blue")
spheres3d(combined_matrix0, radius = 0.001, col = "Red")

# Method 1 

# Use of landmarks to visualise differences using first method
mesh8_warped <-plotRefToTarget(combined_matrix8, combined_matrix0, mesh = mesh8, method = "surface" )

# Compare the two meshes
shade3d(mesh8_warped, col = 2, specular = 1)
shade3d(mesh8, col = 3, specular = 1)

# Visualize within the heatmaps
MD8 <- meshDist(mesh8, mesh8_warped, from=-30, to = 15)

# Method 2 (Less Computationally Demanding)

# Warp the mesh to the average landmark positions
mesh8_warped <- tps3d(mesh8, combined_matrix8, combined_matrix0, threads = 1)

# Compare the two meshes
shade3d(mesh8_warped, col = 2, specular = 1)
shade3d(mesh8, col = 3, specular = 1)

# Visualize within the heatmaps
MD8 <- meshDist(mesh8, mesh8_warped, from = -50, to = 10)



