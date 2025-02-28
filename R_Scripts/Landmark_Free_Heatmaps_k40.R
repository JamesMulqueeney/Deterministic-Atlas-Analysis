# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Kernel 40.0 - Landmark-free Heatmaps 

################################################################################

# Load required libraries 
library(Morpho)
library(abind)
library(rgl)
library(Rvcg)
library(dplyr)
library(geomorph)

#########################################################################################

# Kernel 40.00mm (45 Control Points)

#########################################################################################

# Define input directory
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Outputs/"

################################################################################

# Load in the atlas 
atlas_40.0 <- read.ply(file.path(input_dir, "Binturong_DeterministicAtlas_Kernel_40.ply"))

#########################################################################################

# Apply to Artiodactyla Specimen - Bos domestic

#########################################################################################

# Load Comparison Mesh
mesh1 <- read.ply(file.path(input_dir, "Bos_domestic.ply"))

# Plot warped vs original 
shade3d(mesh1, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.1 <-  meshDist(mesh1, atlas_40.0, from = -150, to = 20)

#########################################################################################

# Apply to Carnivora - Arctictis binturong

#########################################################################################

# Load Comparison Mesh 
mesh2 <- read.ply(file.path(input_dir, "Arctictis_binturong.ply"))

# Plot warped vs original 
shade3d(mesh2, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.2 <-  meshDist(mesh2, atlas_40.0, from = -5, to = 5)

#########################################################################################

# Apply to Cetacea - Schizodelphis morckhoviensis

#########################################################################################

# Load Comparison Mesh  
mesh3 <- read.ply(file.path(input_dir, "Schizodelphis_morckhoviensis.ply"))

# Plot warped vs original 
shade3d(mesh3, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.3 <-  meshDist(mesh3, atlas_40.0, from=-80, to = 20)

#########################################################################################

# Apply to Chiroptera - Centurio senex

#########################################################################################

# Load Comparison Mesh  
mesh4 <- read.ply(file.path(input_dir, "Centurio_senex.ply"))

# Plot warped vs original 
shade3d(mesh4, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.4 <-  meshDist(mesh4, atlas_40.0, from=-40, to = 15)

#########################################################################################

# Apply to Perrisodactyla - Equus caballus

#########################################################################################

# Load Comparison Mesh  
mesh5 <- read.ply(file.path(input_dir, "Equus_caballus.ply"))

# Plot warped vs original 
shade3d(mesh5, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.5 <-  meshDist(mesh5, atlas_40.0, from=-40, to = 15)

#########################################################################################

# Apply to Pilosa - Choloepus hoffmanni

#########################################################################################

# Load Comparison Mesh  
mesh6 <- read.ply(file.path(input_dir, "Choloepus_hoffmanni.ply"))

# Plot warped vs original 
shade3d(mesh6, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.6 <-  meshDist(mesh6, atlas_40.0,  from=-25, to = 15)

#########################################################################################

# Apply to Primates - Cacajao_calvus

#########################################################################################

# Load Comparison Mesh  
mesh7 <- read.ply(file.path(input_dir, "Cacajao_calvus.ply"))

# Plot warped vs original 
shade3d(mesh7, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.7 <-  meshDist(mesh7, atlas_40.0, from=-40, to = 15)

#########################################################################################

# Apply to Rodentia - Cavia aperea

#########################################################################################

# Load Comparison Mesh  
mesh8 <- read.ply(file.path(input_dir, "Cavia_aperea.ply"))

# Plot warped vs original 
shade3d(mesh8, col="Blue", specular=1, alpha = 0.5)
shade3d(atlas_40.0, col="Red", specular=1, alpha = 0.5)

# Plot the warped mesh 
MD_40.8 <-  meshDist(mesh8, atlas_40.0, from=-50, to = 10)
