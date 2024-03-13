# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Produce Heatmaps 

# Author: James M. Mulqueeney 

# Date Last Modified: 06/02/2024 

# Landmark Free Data Only 

# Load in the correct libraries 
library(Morpho)
library(rgl)
library(Rvcg)
library(dplyr)

#########################################################################################
# Use atlas mesh 

# Load in the Original Binturong Specimen 
mesh2 <- read.ply("E:path/to/input/Aligned Meshes/Processed Meshes/ASCII/Cacajao_calvus.ply")

# Scale to the same size as the manual landmark data 
mesh2.scaled <- scalemesh(mesh2, size=0.00044342731, center="none")

#########################################################################################
# Kernel 40.00mm (45 Control Points)

# Load in the atlas 
atlas_40.0 <- read.ply("path/to/input/Poisson Meshes/Kernel 40.0/DeterministicAtlas__EstimatedParameters__Template_cranium.ply")

# Scale to the same size as the manual landmark data 
atlas_40.0.scaled <- scalemesh(atlas_40.0, size=0.00044342731, center="none")

# Plot warped vs original 
shade3d(mesh2.scaled, col=2, specular=2)
shade3d(atlas_40.0.scaled, col=3, specular=2)

# Plot the warped mesh 
MD_40 <-  meshDist(mesh2.scaled, atlas_40.0.scaled, from =-0.02, to=0.01)

#########################################################################################
# Kernel 20.00mm (45 Control Points)

# Load in the atlas 
atlas_20.0 <- read.ply("path/to/input/Poisson Meshes/Kernel 20.0/DeterministicAtlas__EstimatedParameters__Template_cranium.ply")

# Scale to the same size as the manual landmark data 
atlas_20.0.scaled <- scalemesh(atlas_20.0, size=0.00044342731, center="none")

# Plot warped vs original 
shade3d(mesh2.scaled, col=2, specular=2)
shade3d(atlas_20.0.scaled, col=3, specular=2)

# Plot the warped mesh 
MD_20 <-  meshDist(mesh2.scaled, atlas_20.0.scaled, from =-0.02, to=0.01)

#########################################################################################
# Kernel 10.00mm (45 Control Points)

# Load in the atlas 
atlas_10.0 <- read.ply("path/to/input/Kernel 10.0/DeterministicAtlas__EstimatedParameters__Template_cranium.ply")

# Scale to the same size as the manual landmark data 
atlas_10.0.scaled <- scalemesh(atlas_10.0, size=0.00044342731, center="none")

# Plot warped vs original 
shade3d(mesh2.scaled, col=2, specular=2)
shade3d(atlas_10.0.scaled, col=3, specular=2)

# Plot the warped mesh 
MD_10 <-  meshDist(mesh2.scaled, atlas_10.0.scaled, from =-0.02, to=0.01)


