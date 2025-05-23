# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Single Mesh Alignment 

################################################################################

# Load required libraries
library(Morpho)
library(geomorph)
library(abind)
library(rgl)
library(Rvcg)

#########################################################################################

# Read in the landmark data 
data <- read.csv("path/to/input/Mirrored.322.csv", header = TRUE)

# Read in mesh data 
mesh1 <- read.ply("path/to/input/ASCII Mesh Files/Acinonyx_jubatus.ply",  ShowSpecimen= FALSE)
mesh2 <- read.ply("path/to/input/Acratocnus_odontrigonus.ply",  ShowSpecimen= FALSE)
mesh3 <- read.ply("path/to/input/Adapis_magnus.ply",  ShowSpecimen= FALSE)


# Extract x,y,z co-ordinates in 3D
x_values <- as.matrix(data[, grep("\\.X", colnames(data))])
y_values <- as.matrix(data[, grep("\\.Y", colnames(data))])
z_values <- as.matrix(data[, grep("\\.Z", colnames(data))])

num_landmarks <- ncol(x_values)
num_specimens <- nrow(data)

result_matrix <- array(0, dim = c(num_landmarks, 3, num_specimens))

# Fill the result_matrix with x, y, and z values for each specimen
for (i in 1:num_specimens) {
  result_matrix[, 1, i] <- x_values[i, ]
  result_matrix[, 2, i] <- y_values[i, ]
  result_matrix[, 3, i] <- z_values[i, ]
}

# Combine the x,y,z values in a new array 
concatenated_array <- abind(result_matrix, along = 1)

# Extract landmarks for mesh 1
specimen1 <- concatenated_array [,,1]
# Extract landmarks for mesh 2
specimen2 <- concatenated_array [,,2]
# Extract landmarks for mesh 3
specimen3 <- concatenated_array [,,3]

# Plot landmarks with original mesh 

# Mesh 1
spheres3d(specimen1) # Plots the landmarks 
shade3d(mesh1, col=2, specular=1) #Plots the mesh 
# Mesh 2
spheres3d(specimen2) # Plots the landmarks 
shade3d(mesh2, col=3, specular=1) #Plots the mesh 
# Mesh 3
spheres3d(specimen3) # Plots the landmarks 
shade3d(mesh3, col=4, specular=1) #Plots the mesh 

# Mesh alignment

# Plot the first mesh for alignment 
rotmesh1 <-rotmesh.onto(mesh1,specimen1, specimen1, scale=TRUE)
shade3d(rotmesh1$mesh, col=2, specular=1)

# Align the 2nd mesh 
rotmesh2 <-rotmesh.onto(mesh2,specimen2, specimen1, scale=TRUE)
shade3d(rotmesh2$mesh, col=3, specular=1)

# Align the 3rd mesh 
rotmesh3 <-rotmesh.onto(mesh3,specimen3, specimen1, scale=TRUE)
shade3d(rotmesh3$mesh, col=4, specular=1)

# Define the file name for the output mesh
output_file <- "path/to/output/rotmesh3.ply"

# Export the mesh in ASCII format
vcgPlyWrite(rotmesh3$mesh, output_file)
