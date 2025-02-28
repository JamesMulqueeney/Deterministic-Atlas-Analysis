# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Batch Mesh Decimation 

################################################################################

# Load in the correct libraries 
library(Morpho)
library(geomorph)
library(abind)
library(rgl)
library(Rvcg)

################################################################################

# Read in the Mesh files 

################################################################################

# Define the directory containing the .ply files
directory <- "mesh/input/directory"

# Get a list of all .ply files in the directory
ply_files <- list.files(directory, pattern = "\\.ply$", full.names = TRUE)

# Initialize a list to store the mesh objects
mesh_list <- list()

# Extract original file names
original_file_names <- basename(ply_files)

# Loop through the list of .ply files and read them
for (i in seq_along(ply_files)) {
  mesh <- read.ply(ply_files[i], ShowSpecimen= FALSE)
  mesh_list[[i]] <- mesh
}

################################################################################

# Decimate the Mesh Files 

################################################################################
# Initialize a list to store the decimated meshes
decimated_mesh_list <- list()

# Loop through each mesh and apply the decimation
for (i in seq_along(mesh_list)) {
  decimated_mesh_list[[i]] <- vcgQEdecim(mesh_list[[i]], tarface = 50,000)
}

# Check results
print(decimated_mesh_list)
