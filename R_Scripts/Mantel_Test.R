# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Mantel Test Comparisons  

################################################################################

# Load required libraries
library(geomorph)
library(dplyr)
library(vegan)

################################################################################

# Define input directory
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/"

# Read species data
species.data <- read.csv(file.path(input_dir, "Data_A20-Specimen_Details.csv"))

################################################################################

# Manual landmarking Results - Goswami et al. (2022)

################################################################################

# Read in manual landmarking data 
shape.data <- read.csv(file.path(input_dir, "Data_A21-Shape_Data_322.csv"))

# Perform PCA on shape data 
PCA1 <- gm.prcomp(shape.data[,-1])

# Extract summary
PCAsummary<-summary(PCA1)

# Extract the PCA values 
pcscores<-PCA1$x

# Join PCA and Species data together 
pcscores<-cbind(species.data,pcscores) 

################################################################################

# Compare against Aligned-only dataset

################################################################################

# Landmark-free data - Aligned_Only

# Kernel 40.0 (Control points = 45 )
Kernel_40.0 <- read.csv(file.path(input_dir, "Data_A4-Aligned_Only_k40_kpca.csv"))

# Kernel 20.0 (Control points = 270 )
Kernel_20.0 <- read.csv(file.path(input_dir, "Data_A5-Aligned_Only_k20_kpca.csv"))

# Kernel 10.0 (Control points = 1782 )
Kernel_10.0 <- read.csv(file.path(input_dir, "Data_A6-Aligned_Only_k10_kpca.csv"))

################################################################################

# Aligned_Only: Manual vs Kernel 40.0 (Control points = 45 )

################################################################################

# Merge the data into a new dataframe
merged_data1 <- merge(pcscores, Kernel_40.0, by = "Tip_Label")

# Extract relevant components and convert to matrices
pcscores_matrix1 <- as.matrix(merged_data1 %>% select(starts_with("Comp")))
kernel_matrix1 <- as.matrix(merged_data1 %>% select(starts_with("PC")))

# Compute distance matrices
pcscores_dist1 <- dist(pcscores_matrix1)  
kernel_dist1 <- dist(kernel_matrix1)

# Apply Mantel Test
mantel1 <- mantel(pcscores_dist1, kernel_dist1, method = "pearson", permutations = 9999)
print(mantel1)

################################################################################

# Aligned_Only: Manual vs Kernel 20.0 (Control points = 270 )

################################################################################

# Merge the data into a new dataframe 
merged_data2 <- merge(pcscores, Kernel_20.0, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix2 <- as.matrix(merged_data2 %>% select(starts_with("Comp")))
kernel_matrix2 <- as.matrix(merged_data2 %>% select(starts_with("PC")))

# Compute distance matrices
pcscores_dist2 <- dist(pcscores_matrix2)  
kernel_dist2 <- dist(kernel_matrix2)

# Apply Mantel Test
mantel2 <-mantel(pcscores_dist2, kernel_dist2, method = "pearson", permutations = 9999)
print(mantel2)

################################################################################

# Aligned_Only: Manual vs Kernel 10.0 (Control points = 1782 )

################################################################################

# Merge the data into a new dataframe 
merged_data3 <- merge(pcscores, Kernel_10.0, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix3 <- as.matrix(merged_data3 %>% select(starts_with("Comp")))
kernel_matrix3 <- as.matrix(merged_data3 %>% select(starts_with("PC")))

# Compute distance matrices
pcscores_dist3 <- dist(pcscores_matrix3)  
kernel_dist3 <- dist(kernel_matrix3)

# Apply Mantel Test
mantel3 <-mantel(pcscores_dist3, kernel_dist3, method = "pearson", permutations = 9999)
print(mantel3)

################################################################################

# Compare against Poisson mesh dataset

################################################################################

# Landmark-free data - Poisson Meshes 

# Kernel 40.0 (Control points = 45 )
Kernel_40.0_P <- read.csv(file.path(input_dir, "Data_A7-Poisson_k40_kpca.csv"))

# Kernel 20.0 (Control points = 270 )
Kernel_20.0_P <- read.csv(file.path(input_dir, "Data_A8-Poisson_k20_kpca.csv"))

# Kernel 10.0 (Control points = 1782 )
Kernel_10.0_P <- read.csv(file.path(input_dir, "Data_A9-Poisson_k10_kpca.csv"))

################################################################################

# Poisson Mesh: Manual vs Kernel 40.0 (Control points = 45)

################################################################################

# Merge the data into a new dataframe 
merged_data4 <- merge(pcscores, Kernel_40.0_P, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix4 <- as.matrix(merged_data4 %>% select(starts_with("Comp")))
kernel_matrix4 <- as.matrix(merged_data4 %>% select(starts_with("PC")))

# Compute distance matrices
pcscores_dist4 <- dist(pcscores_matrix4)  
kernel_dist4 <- dist(kernel_matrix4)

# Apply Mantel Test
mantel4 <-mantel(pcscores_dist4, kernel_dist4, method = "pearson", permutations = 9999)
print(mantel4)

################################################################################

# Poisson Mesh: Manual vs Kernel 20.0 (Control points = 270)

################################################################################

# Merge the data into a new dataframe 
merged_data5 <- merge(pcscores, Kernel_20.0_P, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix5 <- as.matrix(merged_data5 %>% select(starts_with("Comp")))
kernel_matrix5 <- as.matrix(merged_data5 %>% select(starts_with("PC")))

# Compute distance matrices
pcscores_dist5 <- dist(pcscores_matrix5)  
kernel_dist5 <- dist(kernel_matrix5)

# Apply Mantel Test
mantel5 <-mantel(pcscores_dist5, kernel_dist5, method = "pearson", permutations = 9999)
print(mantel5)

################################################################################

# Poisson Mesh: Manual vs Kernel 10.0 (Control points = 1782)

################################################################################

# Merge the data into a new dataframe 
merged_data6 <- merge(pcscores, Kernel_10.0_P, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix6 <- as.matrix(merged_data6 %>% select(starts_with("Comp")))
kernel_matrix6 <- as.matrix(merged_data6 %>% select(starts_with("PC")))

# Compute distance matrices
pcscores_dist6 <- dist(pcscores_matrix6)  
kernel_dist6 <- dist(kernel_matrix6)

# Apply Mantel Test
mantel6 <-mantel(pcscores_dist6, kernel_dist6, method = "pearson", permutations = 9999)
print(mantel6)

################################################################################

# T-Test Comparison of Aligned-only vs Poisson Results

# Define Mantel statistic values
aligned <- c(0.2799, 0.2647, 0.3098) # Mean =  0.285 SD = 0.0229
poisson <- c(0.6574, 0.6107, 0.5879) # Mean =  0.619 SD = 0.0354

# Measure Means  
signif(mean(aligned), 3)
signif(sd(aligned), 3)

# Measure Standard Deviation 
signif(mean(poisson), 3)
signif(sd(poisson), 3)

# Perform paired t-test
t.test(poisson, aligned, paired = TRUE)

