# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Produce Heatmaps 

# Author: James M. Mulqueeney 

# Date Last Modified: 1/03/2024

# Manual Landmarking Data Only 

# Load in the correct libraries 
library(Morpho)
library(abind)
library(rgl)
library(Rvcg)
library(dplyr)

#########################################################################################

# Load in the Atlas Mesh 
mesh1 <- read.ply("path/to/input/Original Files/ASCII Mesh Files/Decimated/Cacajao_calvus.ply")

#########################################################################################

# Manual Landmarking Heatmaps

# Read in the landmark data 
shape_data <- read.csv("E:/CTData/James Mulqueeney/Mammalian Data/Full Results/Manual Landmarking/Data/shape.data.322.csv", header = TRUE)

# Extract x,y,z co-ordinates in 3D
x_values <- as.matrix(shape_data[, grep("\\.X", colnames(shape_data))])
y_values <- as.matrix(shape_data[, grep("\\.Y", colnames(shape_data))])
z_values <- as.matrix(shape_data[, grep("\\.Z", colnames(shape_data))])

num_landmarks <- ncol(x_values)
num_specimens <- nrow(shape_data)

result_matrix <- array(0, dim = c(num_landmarks, 3, num_specimens))

# Fill the result_matrix with x, y, and z values for each specimen
for (i in 1:num_specimens) {
  result_matrix[, 1, i] <- x_values[i, ]
  result_matrix[, 2, i] <- y_values[i, ]
  result_matrix[, 3, i] <- z_values[i, ]
}

# Combine the x,y,z values in a new array 
shape_concatenated_array <- abind(result_matrix, along = 1)

# Calculate the average for each column across all matrix slices
shape_average_matrix_slice <- apply(shape_concatenated_array, MARGIN = c(1, 2), FUN = mean)

# Compare the landmark positions 

# Plot the average_matrix_slice
spheres3d(shape_average_matrix_slice, radius=0.001, col="Red")

# Plot against the atlas 
spheres3d(shape_concatenated_array [,,36], radius=0.001) #Arctictis is No.22 in dataset

# Visualise the differences between landmarks 
deformGrid3d(shape_concatenated_array [,,36],shape_average_matrix_slice,ngrid = 5)

#########################################################################################

# Repeat for the mirrored data (minus the mirroring)

# Read in original positions 
mirrored_data <- read.csv("path/to/input/Mirrored.322_Reduced.csv", header = TRUE)

# Extract x,y,z co-ordinates in 3D
x_values <- as.matrix(mirrored_data[, grep("\\.X", colnames(mirrored_data))])
y_values <- as.matrix(mirrored_data[, grep("\\.Y", colnames(mirrored_data))])
z_values <- as.matrix(mirrored_data[, grep("\\.Z", colnames(mirrored_data))])

num_landmarks <- ncol(x_values)
num_specimens <- nrow(mirrored_data)

result_matrix <- array(0, dim = c(num_landmarks, 3, num_specimens))

# Fill the result_matrix with x, y, and z values for each specimen
for (i in 1:num_specimens) {
  result_matrix[, 1, i] <- x_values[i, ]
  result_matrix[, 2, i] <- y_values[i, ]
  result_matrix[, 3, i] <- z_values[i, ]
}

# Combine the x,y,z values in a new array 
mirrored_concatenated_array <- abind(result_matrix, along = 1)

# Calculate the average for each column across all matrix slices
mirrored_average_matrix_slice <- apply(mirrored_concatenated_array, MARGIN = c(1, 2), FUN = mean)

# Compare the landmark positions 

# Plot the average_matrix_slice
spheres3d(mirrored_average_matrix_slice, col="Red")

# Plot against the atlas 
spheres3d(mirrored_concatenated_array [,,36]) #Arctictis is No.22 in dataset

# Deform the mesh into the correct size 
mesh1.corrected <- tps3d(mesh1, mirrored_concatenated_array [,,36], shape_concatenated_array [,,36])
vcgPlyWrite(mesh1.corrected, "deformed_mesh")

# Check with the landmarks that it has been correctly deformed 
shade3d(mesh1.corrected, col=3, specular=1)
spheres3d(shape_concatenated_array [,,36], radius=0.001)

#########################################################################################

# Apply to specimen 

# Warp mesh to the the average landmark positions 
mesh1.warped <- tps3d(mesh1.corrected,shape_concatenated_array [,,36],shape_average_matrix_slice,threads=1)

# Compare the two meshes with each other 
shade3d(mesh1.warped, col=2, specular=1)
spheres3d(shape_average_matrix_slice, radius=0.001, col=2)

shade3d(mesh1.corrected, col=3, specular=1)
spheres3d(shape_concatenated_array [,,22], radius=0.001, col=3)

# Now use to visualise within the heatmaps 
MD1 <-  meshDist(mesh1.corrected, mesh1.warped, from =-0.02, to=0.01)

