# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Create Evolutionary Rates & Disparity Data using Momentas 

################################################################################

# Read in libraries 
library(geomorph)
library(ape)

################################################################################

# Kernel 40.0 (Control points = 45) 

################################################################################

# Read the first line of the file to get the dimensions
k40_first_line <- readLines("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Large Results/Final Poisson Analysis/Arctictis_binturong_atlas/Kernel 40.0/output/DeterministicAtlas__EstimatedParameters__Momenta.txt", n=1)

# Split the first line into individual elements
k40_dimensions <- as.numeric(strsplit(k40_first_line, " ")[[1]])

# Assign dimensions to variables
k40_num_specimens <- k40_dimensions[1]
k40_lines_per_specimen <- k40_dimensions[2]
k40_num_values_per_line <- k40_dimensions[3]

# Read the data file skipping the first line
k40_data <- scan("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Large Results/Final Poisson Analysis/Arctictis_binturong_atlas/Kernel 40.0/output/DeterministicAtlas__EstimatedParameters__Momenta.txt", skip = 1)

# Reshape the data into a 3D array
k40_data_array <- array(k40_data, dim = c(k40_num_values_per_line, k40_lines_per_specimen, k40_num_specimens))

# Initialize a list to store each specimen's data
k40_specimens <- vector("list", k40_num_specimens)

# Loop through each specimen and extract its data into a 2D array
for (i in 1:k40_num_specimens) {
  k40_specimens[[i]] <- t(k40_data_array[,,i])
}

# Now you can access each specimen separately
# For example, to access the first specimen:
k40_first_specimen <- k40_specimens[[1]]

# Convert the list of specimen arrays into a 3D array
k40_3d_array <- array(unlist(k40_specimens), dim = c(k40_lines_per_specimen, k40_num_values_per_line, k40_num_specimens))

################################################################################

# Read in the CSV file
Kernel_40.0 <- read.csv("C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/Data_A7-Poisson_k40_kpca.csv")

# Assign the tip labels to the matrix 
dimnames(k40_3d_array) <- list(NULL,  NULL, as.character(Kernel_40.0$Tip_Label))

# Construct the geomorph data frame with appropriate checks
gdf_k40 <- geomorph.data.frame(
  coords = k40_3d_array,
  specimen_id = as.factor(Kernel_40.0$Tip_Label),  # Ensure factor type for grouping variables
  diet = as.factor(Kernel_40.0$Diet), 
  loc = as.factor(Kernel_40.0$Locomotion), 
  dev = as.factor(Kernel_40.0$Development), 
  clade = as.factor(Kernel_40.0$Order), 
  superorder = as.factor(Kernel_40.0$Superorder)
)

################################################################################

# Phylogenetic Signal 

################################################################################

# Read in the Dated Tree File 
phy <- read.tree("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Deterministic-Atlas-Analysis-main (2)/Deterministic-Atlas-Analysis-main/R_Scripts/Adapted Code/trees90_95subset.tre")

# Extract the first tree 
single_phy <- phy[[2]] 

# Run the Phylogenetic Signal Analysis 
physig_k40 <- physignal(A = gdf_k40$coords, phy = single_phy, iter = 99)
print(physig_k40)

################################################################################

# Disparity 

################################################################################

# Disparity by Diet 
disparity_diet_k40 <- morphol.disparity(f1 = coords ~ diet, groups = gdf_k40$diet, data = gdf_k40)
print(disparity_diet_k40)

# Disparity by Locomotion 
disparity_loc_k40 <- morphol.disparity(f1 = coords ~ loc, groups = gdf_k40$loc, data = gdf_k40)
print(disparity_loc_k40)

# Disparity by Development 
disparity_dev_k40 <- morphol.disparity(f1 = coords ~ dev, groups = gdf_k40$dev, data = gdf_k40)
print(disparity_dev_k40)

# Disparity by Order (clade)
disparity_clade_k40 <- morphol.disparity(f1 = coords ~ clade, groups = gdf_k40$clade, data = gdf_k40)
print(disparity_clade_k40)

# Disparity by Superorder 
disparity_super_k40 <- morphol.disparity(f1 = coords ~ superorder, groups = gdf_k40$superorder, data = gdf_k40)
print(disparity_super_k40)

################################################################################

# Evolutionary rates

################################################################################

# Rates by Diet 
gp.diet <- gdf_k40$diet
names(gp.diet) <- gdf_k40$specimen_id

dietrates_k40 <- compare.evol.rates(A = gdf_k40$coords, gp = gp.diet, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(dietrates_k40)

# Rates by Locomotion 
gp.loc <- gdf_k40$loc
names(gp.loc) <- gdf_k40$specimen_id

locrates_k40 <- compare.evol.rates(A = gdf_k40$coords, gp = gp.loc, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(locrates_k40)

# Rates by Development 
gp.dev <- gdf_k40$dev
names(gp.dev) <- gdf_k40$specimen_id

devorates_k40 <- compare.evol.rates(A = gdf_k40$coords, gp = gp.dev, phy = single_phy, iter = 100, method = "simulation", print.progress = TRUE)
print(devorates_k40)

# Rates by Order (clade)
gp.class <- gdf_k40$clade
names(gp.class) <- gdf_k40$specimen_id

classrates_k40 <- compare.evol.rates(A = gdf_k40$coords, gp = gp.class, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(classrates_k40)

# Rates by Superorder 
gp.super <- gdf_k40$superorder
names(gp.super) <- gdf_k40$specimen_id

superorderrates_k40 <- compare.evol.rates(A = gdf_k40$coords, gp = gp.super, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(superorderrates_k40)

################################################################################

# # Kernel 20.0 (Control points = 270)

################################################################################

# Read the first line of the file to get the dimensions
k20_first_line <- readLines("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Large Results/Final Poisson Analysis/Arctictis_binturong_atlas/Kernel 20.0/output/DeterministicAtlas__EstimatedParameters__Momenta.txt", n=1)

# Split the first line into individual elements
k20_dimensions <- as.numeric(strsplit(k20_first_line, " ")[[1]])

# Assign dimensions to variables
k20_num_specimens <- k20_dimensions[1]
k20_lines_per_specimen <- k20_dimensions[2]
k20_num_values_per_line <- k20_dimensions[3]

# Read the data file skipping the first line
k20_data <- scan("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Large Results/Final Poisson Analysis/Arctictis_binturong_atlas/Kernel 20.0/output/DeterministicAtlas__EstimatedParameters__Momenta.txt", skip = 1)

# Reshape the data into a 3D array
k20_data_array <- array(k20_data, dim = c(k20_num_values_per_line, k20_lines_per_specimen, k20_num_specimens))

# Initialize a list to store each specimen's data
k20_specimens <- vector("list", k20_num_specimens)

# Loop through each specimen and extract its data into a 2D array
for (i in 1:k20_num_specimens) {
  k20_specimens[[i]] <- t(k20_data_array[,,i])
}

# Now you can access each specimen separately
# For example, to access the first specimen:
k20_first_specimen <- k20_specimens[[1]]

# Convert the list of specimen arrays into a 3D array
k20_3d_array <- array(unlist(k20_specimens), dim = c(k20_lines_per_specimen, k20_num_values_per_line, k20_num_specimens))

################################################################################

# Read in the CSV file
Kernel_20.0 <- read.csv("C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/Data_A8-Poisson_k20_kpca.csv")

# Assign the tip labels to the matrix 
dimnames(k20_3d_array) <- list(NULL,  NULL, as.character(Kernel_20.0$Tip_Label))

# Construct the geomorph data frame with appropriate checks
gdf_k20 <- geomorph.data.frame(
  coords = k20_3d_array,
  specimen_id = as.factor(Kernel_20.0$Tip_Label),  # Ensure factor type for grouping variables
  diet = as.factor(Kernel_20.0$Diet), 
  loc = as.factor(Kernel_20.0$Locomotion), 
  dev = as.factor(Kernel_20.0$Development), 
  clade = as.factor(Kernel_20.0$Order), 
  superorder = as.factor(Kernel_20.0$Superorder)
)

# Run the Phylogenetic Signal Analysis 
physig_k20 <- physignal(A = gdf_k20$coords, phy = single_phy, iter = 99)
print(physig_k20)

################################################################################

# Disparity 

################################################################################

# Disparity by Diet 
disparity_diet_k20 <- morphol.disparity(f1 = coords ~ diet, groups = gdf_k20$diet, data = gdf_k20)
print(disparity_diet_k20)

# Disparity by Locomotion 
disparity_loc_k20 <- morphol.disparity(f1 = coords ~ loc, groups = gdf_k20$loc, data = gdf_k20)
print(disparity_loc_k20)

# Disparity by Development 
disparity_dev_k20 <- morphol.disparity(f1 = coords ~ dev, groups = gdf_k20$dev, data = gdf_k20)
print(disparity_dev_k20)

# Disparity by Order (clade)
disparity_clade_k20 <- morphol.disparity(f1 = coords ~ clade, groups = gdf_k20$clade, data = gdf_k20)
print(disparity_clade_k20)

# Disparity by Superorder 
disparity_super_k20 <- morphol.disparity(f1 = coords ~ superorder, groups = gdf_k20$superorder, data = gdf_k20)
print(disparity_super_k20)

################################################################################

# Evolutionary rates

################################################################################

# Rates by Diet 
gp.diet <- gdf_k20$diet
names(gp.diet) <- gdf_k20$specimen_id

dietrates_k20 <- compare.evol.rates(A = gdf_k20$coords, gp = gp.diet, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(dietrates_k20)

# Rates by Locomotion 
gp.loc <- gdf_k20$loc
names(gp.loc) <- gdf_k20$specimen_id

locrates_k20 <- compare.evol.rates(A = gdf_k20$coords, gp = gp.loc, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(locrates_k20)

# Rates by Development 
gp.dev <- gdf_k20$dev
names(gp.dev) <- gdf_k20$specimen_id

devorates_k20 <- compare.evol.rates(A = gdf_k20$coords, gp = gp.dev, phy = single_phy, iter = 100, method = "simulation", print.progress = TRUE)
print(devorates_k20)

# Rates by Order (clade)
gp.class <- gdf_k20$clade
names(gp.class) <- gdf_k20$specimen_id

classrates_k20 <- compare.evol.rates(A = gdf_k20$coords, gp = gp.class, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(classrates_k20)

# Rates by Superorder 
gp.super <- gdf_k20$superorder
names(gp.super) <- gdf_k20$specimen_id

superorderrates_k20 <- compare.evol.rates(A = gdf_k20$coords, gp = gp.super, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(superorderrates_k20)

################################################################################

# Kernel 10.0 (Control points = 1782)

################################################################################

# Read the first line of the file to get the dimensions
k10_first_line <- readLines("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Large Results/Final Poisson Analysis/Arctictis_binturong_atlas/Kernel 10.0/output/DeterministicAtlas__EstimatedParameters__Momenta.txt", n=1)

# Split the first line into individual elements
k10_dimensions <- as.numeric(strsplit(k10_first_line, " ")[[1]])

# Assign dimensions to variables
k10_num_specimens <- k10_dimensions[1]
k10_lines_per_specimen <- k10_dimensions[2]
k10_num_values_per_line <- k10_dimensions[3]

# Read the data file skipping the first line
k10_data <- scan("E:/James Mulqueeney/Paper 2- Mammal Shape/Final Results and Code/Large Results/Final Poisson Analysis/Arctictis_binturong_atlas/Kernel 10.0/output/DeterministicAtlas__EstimatedParameters__Momenta.txt", skip = 1)

# Reshape the data into a 3D array
k10_data_array <- array(k10_data, dim = c(k10_num_values_per_line, k10_lines_per_specimen, k10_num_specimens))

# Initialize a list to store each specimen's data
k10_specimens <- vector("list", k10_num_specimens)

# Loop through each specimen and extract its data into a 2D array
for (i in 1:k10_num_specimens) {
  k10_specimens[[i]] <- t(k10_data_array[,,i])
}

# Now you can access each specimen separately
# For example, to access the first specimen:
k10_first_specimen <- k10_specimens[[1]]

# Convert the list of specimen arrays into a 3D array
k10_3d_array <- array(unlist(k10_specimens), dim = c(k10_lines_per_specimen, k10_num_values_per_line, k10_num_specimens))

################################################################################

# Read in the CSV file
Kernel_10.0 <- read.csv("C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/Data_A9-Poisson_k10_kpca.csv")

# Assign the tip labels to the matrix 
dimnames(k10_3d_array) <- list(NULL,  NULL, as.character(Kernel_10.0$Tip_Label))

# Construct the geomorph data frame with appropriate checks
gdf_k10 <- geomorph.data.frame(
  coords = k10_3d_array,
  specimen_id = as.factor(Kernel_10.0$Tip_Label),  # Ensure factor type for grouping variables
  diet = as.factor(Kernel_10.0$Diet), 
  loc = as.factor(Kernel_10.0$Locomotion), 
  dev = as.factor(Kernel_10.0$Development), 
  clade = as.factor(Kernel_10.0$Order), 
  superorder = as.factor(Kernel_10.0$Superorder)
)

# Run the Phylogenetic Signal Analysis 
physig_k10 <- physignal(A = gdf_k10$coords, phy = single_phy, iter = 99)
print(physig_k10)

################################################################################

# Disparity 

################################################################################

# Disparity by Diet 
disparity_diet_k10 <- morphol.disparity(f1 = coords ~ diet, groups = gdf_k20$diet, data = gdf_k10)
print(disparity_diet_k10)

# Disparity by Locomotion 
disparity_loc_k10 <- morphol.disparity(f1 = coords ~ loc, groups = gdf_k10$loc, data = gdf_k10)
print(disparity_loc_k10)

# Disparity by Development 
disparity_dev_k10 <- morphol.disparity(f1 = coords ~ dev, groups = gdf_k10$dev, data = gdf_k10)
print(disparity_dev_k10)

# Disparity by Order (clade)
disparity_clade_k10 <- morphol.disparity(f1 = coords ~ clade, groups = gdf_k20$clade, data = gdf_k10)
print(disparity_clade_k10)

# Disparity by Superorder 
disparity_super_k10 <- morphol.disparity(f1 = coords ~ superorder, groups = gdf_k40$superorder, data = gdf_k10)
print(disparity_super_k10)

################################################################################

# Evolutionary rates

################################################################################

# Rates by Diet 
gp.diet <- gdf_k10$diet
names(gp.diet) <- gdf_k10$specimen_id

dietrates_k10  <- compare.evol.rates(A = gdf_k10$coords, gp = gp.diet, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(dietrates_k10)

# Rates by Locomotion 
gp.loc <- gdf_k10$loc
names(gp.loc) <- gdf_k10$specimen_id

locrates_k10  <- compare.evol.rates(A = gdf_k10$coords, gp = gp.loc, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(locrates_k10)

# Rates by Development 
gp.dev <- gdf_k10$dev
names(gp.dev) <- gdf_k10$specimen_id

devorates_k10  <- compare.evol.rates(A = gdf_k10$coords, gp = gp.dev, phy = single_phy, iter = 100, method = "simulation", print.progress = TRUE)
print(devorates_k10)

# Rates by Order (clade)
gp.class <- gdf_k10$clade
names(gp.class) <- gdf_k10$specimen_id

classrates_k10  <- compare.evol.rates(A = gdf_k10$coords, gp = gp.class, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(classrates_k10)

# Rates by Superorder 
gp.super <- gdf_k10$superorder
names(gp.super) <- gdf_k10$specimen_id

superorderrates_k10 <- compare.evol.rates(A = gdf_k10$coords, gp = gp.super, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(superorderrates_k10)





