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

# Example matrices (replace these with your actual data)
pcscores_matrix1 <- as.matrix(merged_data1 %>% select(starts_with("Comp")))
kernel_matrix1 <- as.matrix(merged_data1 %>% select(starts_with("PC")))

# Perform Procrustes analysis
procrustes_result1 <- procrustes(pcscores_matrix1, kernel_matrix1)

# Print the Procrustes analysis results
print(procrustes_result1)
summary(procrustes_result1)

# Perform permutation test (PROTEST) to assess significance
protest_result1 <- protest(pcscores_matrix1, kernel_matrix1, permutations = 9999)
print(protest_result1)

# Summary of the protest result
summary(protest_result1)

#########################################################################################

# Aligned_Only: Manual vs Kernel 20.0 (Control points = 270 )

#########################################################################################

# Merge the data into a new dataframe 
merged_data2 <- merge(pcscores, Kernel_20.0, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix2 <- as.matrix(merged_data2 %>% select(starts_with("Comp")))
kernel_matrix2 <- as.matrix(merged_data2 %>% select(starts_with("PC")))

# Perform Procrustes analysis
procrustes_result2 <- procrustes(pcscores_matrix2, kernel_matrix2)

# Print the Procrustes analysis results
print(procrustes_result2)
summary(procrustes_result2)

# Perform permutation test (PROTEST) to assess significance
protest_result2 <- protest(pcscores_matrix2, kernel_matrix2, permutations = 9999)
print(protest_result2)

# Summary of the protest result
summary(protest_result2)

#########################################################################################

# Aligned_Only: Manual vs Kernel 10.0 (Control points = 1782 )

#########################################################################################

# Merge the data into a new dataframe 
merged_data3 <- merge(pcscores, Kernel_10.0, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix3 <- as.matrix(merged_data3 %>% select(starts_with("Comp")))
kernel_matrix3 <- as.matrix(merged_data3 %>% select(starts_with("PC")))

# Perform Procrustes analysis
procrustes_result3 <- procrustes(pcscores_matrix3, kernel_matrix3)

# Print the Procrustes analysis results
print(procrustes_result3)
summary(procrustes_result3)

# Perform permutation test (PROTEST) to assess significance
protest_result3 <- protest(pcscores_matrix3, kernel_matrix3, permutations = 9999)
print(protest_result3)

# Summary of the protest result
summary(protest_result3)

#########################################################################################

# Compare against Poisson mesh dataset

#########################################################################################

# Landmark-free data - Poisson Meshes 

# Kernel 40.0 (Control points = 45 )
Kernel_40.0_P <- read.csv(file.path(input_dir, "Data_A7-Poisson_k40_kpca.csv"))

# Kernel 20.0 (Control points = 270 )
Kernel_20.0_P <- read.csv(file.path(input_dir, "Data_A8-Poisson_k20_kpca.csv"))

# Kernel 10.0 (Control points = 1782 )
Kernel_10.0_P <- read.csv(file.path(input_dir, "Data_A9-Poisson_k10_kpca.csv"))


#########################################################################################

# Poisson Mesh: Manual vs Kernel 40.0 (Control points = 45)

#########################################################################################

# Merge the data into a new dataframe 
merged_data4 <- merge(pcscores, Kernel_40.0_P, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix4 <- as.matrix(merged_data4 %>% select(starts_with("Comp")))
kernel_matrix4 <- as.matrix(merged_data4 %>% select(starts_with("PC")))

# Perform Procrustes analysis
procrustes_result4 <- procrustes(pcscores_matrix4, kernel_matrix4)

# Print the Procrustes analysis results
print(procrustes_result4)
summary(procrustes_result4)

# Perform permutation test (PROTEST) to assess significance
protest_result4 <- protest(pcscores_matrix4, kernel_matrix4, permutations = 9999)
print(protest_result4)

# Summary of the protest result
summary(protest_result4)

#########################################################################################

# Poisson Mesh: Manual vs Kernel 20.0 (Control points = 270)

#########################################################################################

# Merge the data into a new dataframe 
merged_data5 <- merge(pcscores, Kernel_20.0_P, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix5 <- as.matrix(merged_data5 %>% select(starts_with("Comp")))
kernel_matrix5 <- as.matrix(merged_data5 %>% select(starts_with("PC")))

# Perform Procrustes analysis
procrustes_result5 <- procrustes(pcscores_matrix5, kernel_matrix5)

# Print the Procrustes analysis results
print(procrustes_result5)
summary(procrustes_result5)

# Perform permutation test (PROTEST) to assess significance
protest_result5 <- protest(pcscores_matrix5, kernel_matrix5, permutations = 9999)
print(protest_result5)

# Summary of the protest result
summary(protest_result5)

#########################################################################################

# Poisson Mesh: Manual vs Kernel 10.0 (Control points = 1782)

#########################################################################################

# Merge the data into a new dataframe 
merged_data6 <- merge(pcscores, Kernel_10.0_P, by = "Tip_Label")

# Example matrices (replace these with your actual data)
pcscores_matrix6 <- as.matrix(merged_data6 %>% select(starts_with("Comp")))
kernel_matrix6 <- as.matrix(merged_data6 %>% select(starts_with("PC")))

# Perform Procrustes analysis
procrustes_result6 <- procrustes(pcscores_matrix6, kernel_matrix6)

# Print the Procrustes analysis results
print(procrustes_result6)
summary(procrustes_result6)

# Perform permutation test (PROTEST) to assess significance
protest_result6 <- protest(pcscores_matrix6, kernel_matrix6, permutations = 9999)
print(protest_result6)

# Summary of the protest result
summary(protest_result6)

########################################################################################

# T-Test Comparison of Aligned-only vs Poisson Results

########################################################################################

# Procrustes sum of Squares 

########################################################################################

# Define Protest statistic values
aligned <- c(0.7619932, 0.6659323, 0.6079318 ) # Mean =  0.679 SD = 0.0778
poisson <- c(0.5030495, 0.4178207, 0.429891 ) # Mean =  0.450 SD = 0.0461

# Measure Means  
signif(mean(aligned), 3)
signif(sd(aligned), 3)

# Measure Standard Deviation 
signif(mean(poisson), 3)
signif(sd(poisson), 3)

# Perform paired t-test
t.test(poisson, aligned, paired = TRUE)


########################################################################################

# Procrustes mean squared error

########################################################################################

# Define Protest statistic values
aligned <- c(0.04864605, 0.04547651, 0.04345097) # Mean =  0.0459 SD = 0.00262
poisson <- c(0.0395255, 0.03602193 , 0.03653855 ) # Mean =  0.0374 SD = 0.00189

# Measure Means  
signif(mean(aligned), 3)
signif(sd(aligned), 3)

# Measure Standard Deviation 
signif(mean(poisson), 3)
signif(sd(poisson), 3)

# Perform paired t-test
t.test(poisson, aligned, paired = TRUE)




