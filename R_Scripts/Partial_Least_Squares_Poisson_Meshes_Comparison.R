# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Partial Least Squares Comparison 

# Author: James M. Mulqueeney 

# Date Last Modified: 13/03/2024

# Poisson Meshes

# Load in the correct libraries 
library(Morpho)
library(abind)
library(rgl)
library(Rvcg)
library(dplyr)
library(geomorph)

#########################################################################################
# Partial Least Squared Analysis

# Manual comparison to all other results 

# Can test variation by PC axes and also across different groups? 

#########################################################################################

# Read in the Data 

# Manual landmarking data 

# Shape data (Landmark Data)
shape.data <- read.csv("path/to/input/Data_S3-Shape_Data_322.csv")

# Species data (Details of Taxonomy etc.)
species.data <- read.csv("path/to/input/Data_S1-Specimen_Details.csv")

# Perform PCA on shape data 
PCA1 <- gm.prcomp(shape.data[,-1])

# Extract summary and write to a .txt file 
PCAsummary<-summary(PCA1)
write.table(PCAsummary$PC.summary, file = "pcasummary.txt")

# Extract the PCA values 
pcscores<-PCA1$x

# Join PCA and Species data together 
pcscores<-cbind(species.data,pcscores) 

#########################################################################################

# Landmark-free data 

# Kernel 40.0 (Control points = 45 )
Kernel_40.0 <- read.csv("path/to/input/Data_S8-Possion_k40_kpca.csv")

# Kernel 20.0 (Control points = 270 )
Kernel_20.0 <- read.csv("path/to/input/Data_S9-Poisson_k20_kpca.csv")

# Kernel 10.0 (Control points = 1782 )
Kernel_10.0 <- read.csv("path/to/input/Data_S10-Poisson_k10_kpca.csv")

#########################################################################################

# For making plots

par(mfrow = c(1, 3)) 

# Manual vs Kernel 40.0 

# Merge the data into a new dataframe 
merged_data1 <- merge(pcscores, Kernel_40.0, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result1 <- two.b.pls(merged_data1[, comps], merged_data1[, pcs])
summary(pls_result1)

# Visualize the results
PLS1<- plot(pls_result1, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

#########################################################################################

# Manual vs Kernel 20.0 

# Merge the data into a new dataframe 
merged_data2 <- merge(pcscores, Kernel_20.0, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result2 <- two.b.pls(merged_data2[, comps], merged_data2[, pcs])
summary(pls_result2)

# Visualize the results
PLS2<- plot(pls_result2, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

#########################################################################################

# Manual vs Kernel 10.0 

# Merge the data into a new dataframe 
merged_data3 <- merge(pcscores, Kernel_10.0, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result3 <- two.b.pls(merged_data3[, comps], merged_data3[, pcs])
summary(pls_result3)

# Visualize the results
PLS3<- plot(pls_result3, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

#########################################################################################

# For all total 
par(mfrow = c(2, 3)) 

# Combined plots 
par(mfrow = c(1, 1))
