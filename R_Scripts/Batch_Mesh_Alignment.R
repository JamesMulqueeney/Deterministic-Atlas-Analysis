# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Alignment of Meshes 

# Author: James M. Mulqueeney 

# Date Last Modified: 13/03/2024

# Single Mesh Alignment 

# Load in the correct libraries 
library(Morpho)
library(geomorph)
library(abind)
library(rgl)
library(Rvcg)

#########################################################################################

# Define the directory containing the .ply files
directory <- "path/to/input/ASCII Mesh Files/"

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

#########################################################################################

# Assuming concatenated_array is a 3D array where the third dimension represents specimens

# Number of specimens
num_specimens <- dim(concatenated_array)[3]

# Initialize a list to store landmarks for each specimen
specimen_list <- vector("list", length = num_specimens)

# Loop through specimens and extract landmarks
for (i in 1:num_specimens) {
  specimen_list[[i]] <- concatenated_array[,,i]
}

# Define a list to store aligned meshes
aligned_meshes <- vector("list", length = length(mesh_list))

# Align each mesh
for (i in seq_along(mesh_list)) {
  aligned_mesh <- rotmesh.onto(mesh_list[[i]], specimen_list[[i]], specimen1, scale = TRUE)
  aligned_meshes[[i]] <- aligned_mesh$mesh
}

# Define a directory to save the aligned meshes
output_directory <- "path/to/input//Aligned Mesh Files/ASCII"

# Set the working directory to the output directory
setwd(output_directory)

#########################################################################################

# Save as non-ascii (preferred) format - binary  

# Loop through aligned meshes, save them by their original file names in binary format
for (i in seq_along(mesh_list)) {
  original_name <- gsub("\\.ply$", "", original_file_names[i])  # remove the .ply extension
  output_file <- paste0(original_name, ".ply")
  # Export the mesh in binary format (choose either BE or LE)
  vcgPlyWrite(aligned_meshes[[i]], output_file, format = "PLY_BINARY_LE")  # Use PLY_BINARY_BE for big-endian
}

#########################################################################################

# Save as ascii format 

# Loop through aligned meshes, save them by their original file names in ASCII format
for (i in seq_along(mesh_list)) {
  original_name <- gsub("\\.ply$", "", original_file_names[i])  # remove the .ply extension
  output_file <- paste0(original_name, ".ply")
  # Export the mesh in ASCII format
  vcgPlyWrite(aligned_meshes[[i]], output_file)
}

