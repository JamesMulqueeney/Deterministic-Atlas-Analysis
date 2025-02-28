# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Create Evolutionary Rates & Disparity Data

################################################################################

# Load required libraries
library(geomorph)
library(ape)

################################################################################

# Manual Landmarking Results 

################################################################################

# Define input directory
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/"

# Read in the CSV file
data <- read.csv(file.path(input_dir, "Data_A24-Combined_Data_322.csv"), stringsAsFactors = FALSE)

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

################################################################################

# Phylogenetic Signal 

################################################################################

# Read in the Dated Tree File 
phy <- read.tree(file.path(input_dir, "trees90_95subset.tre"))

# Extract the first tree 
single_phy <- phy[[2]] 

# Run the Phylogenetic Signal Analysis 
physig <- physignal(A = gdf$coords, phy = single_phy, iter = 99)
print(physig)

################################################################################

# Disparity 

################################################################################

# Disparity by Diet 
disparity_diet <- morphol.disparity(f1 = coords ~ diet, groups = gdf$diet, data = gdf)
print(disparity_diet)

# Disparity by Locomotion 
disparity_loc <- morphol.disparity(f1 = coords ~ loc, groups = gdf$loc, data = gdf)
print(disparity_loc)

# Disparity by Development 
disparity_dev <- morphol.disparity(f1 = coords ~ dev, groups = gdf$dev, data = gdf)
print(disparity_dev)

# Disparity by Order (clade)
disparity_clade <- morphol.disparity(f1 = coords ~ clade, groups = gdf$clade, data = gdf)
print(disparity_clade)

# Disparity by Superorder 
disparity_super <- morphol.disparity(f1 = coords ~ superorder, groups = gdf$superorder, data = gdf)
print(disparity_super)

################################################################################

# Evolutionary rates

################################################################################

# Rates by Diet 
gp.diet <- gdf$diet
names(gp.diet) <- gdf$specimen_id

dietrates <- compare.evol.rates(A =gdf$coords, gp = gp.diet, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(dietrates)

# Rates by Locomotion 
gp.loc <- gdf$loc
names(gp.loc) <- gdf$specimen_id

locrates <- compare.evol.rates(A =gdf$coords, gp = gp.loc, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(locrates)

# Rates by Development 
gp.dev <- gdf$dev
names(gp.dev) <- gdf$specimen_id

devorates <- compare.evol.rates(A = gdf$coords, gp = gp.dev, phy = single_phy, iter = 100, method = "simulation", print.progress = TRUE)
print(devorates)

# Rates by Order (clade)
gp.class <- gdf$clade
names(gp.class) <- gdf$specimen_id

classrates <- compare.evol.rates(A = gdf$coords, gp = gp.class, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(classrates)

# Rates by Superorder 
gp.super <- gdf$superorder
names(gp.super) <- gdf$specimen_id

superorderrates<- compare.evol.rates(A = gdf$coords, gp = gp.super, phy = single_phy, iter=100, method = "simulation", print.progress = TRUE)
print(superorderrates)

################################################################################

# Landmark-free Results 

################################################################################

# Kernel 40.0 (Control points = 45)

################################################################################

# Kernel 40.0 (Control points = 45)
Kernel_40.0 <- read.csv(file.path(input_dir, "Data_A7-Poisson_k40_kpca.csv"))

# Subset the data for shape only 
Kernel_40.0_sub <- Kernel_40.0[19:339]
rownames(Kernel_40.0_sub) <- Kernel_40.0$Tip_Label

# Make into a matrix
Kernel_40.0_matrix <- as.matrix(Kernel_40.0_sub)

# Construct the geomorph data frame with appropriate checks
gdf_k40 <- geomorph.data.frame(
  coords = Kernel_40.0_matrix,
  specimen_id = as.factor(Kernel_40.0$Tip_Label),  # Ensure factor type for grouping variables
  diet = as.factor(Kernel_40.0$Diet), 
  loc = as.factor(Kernel_40.0$Locomotion), 
  dev = as.factor(Kernel_40.0$Development), 
  clade = as.factor(Kernel_40.0$Order), 
  superorder = as.factor(Kernel_40.0$Superorder)
)

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

# Kernel 20.0 (Control points = 270)

################################################################################

# Kernel 20.0 (Control points = 270)
Kernel_20.0 <- read.csv(file.path(input_dir, "Data_A8-Poisson_k20_kpca.csv"))

# Subset the data for shape only 
Kernel_20.0_sub <- Kernel_20.0[19:339]
rownames(Kernel_20.0_sub) <- Kernel_20.0$Tip_Label

# Make into a matrix
Kernel_20.0_matrix <- as.matrix(Kernel_20.0_sub)

# Construct the geomorph data frame with appropriate checks
gdf_k20 <- geomorph.data.frame(
  coords = Kernel_20.0_matrix,
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

# Kernel 10.0 (Control points = 1782)
Kernel_10.0 <- read.csv(file.path(input_dir, "Data_A9-Poisson_k10_kpca.csv"))

# Subset the data for shape only 
Kernel_10.0_sub <- Kernel_10.0[19:339]
rownames(Kernel_10.0_sub) <- Kernel_10.0$Tip_Label

# Make into a matrix
Kernel_10.0_matrix <- as.matrix(Kernel_10.0_sub)

# Construct the geomorph data frame with appropriate checks
gdf_k10 <- geomorph.data.frame(
  coords = Kernel_10.0_matrix,
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





