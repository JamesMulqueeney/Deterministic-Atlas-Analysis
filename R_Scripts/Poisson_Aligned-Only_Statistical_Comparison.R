# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Partial Least Squares Comparison 

# Author: James M. Mulqueeney 

# Date Last Modified: 06/02/2024

# Load in the correct libraries 
library(reshape2)
library(car)

# ANOVA Test 

# Read in PLS Results 
df<- read.csv("path/to/input/Data_S18-Partial_Least_Sqaures.csv")

# Rename columns to remove spaces
colnames(df) <- c("Mesh_Type", "Control_Points_45", "Control_Points_270", "Control_Points_1782")

# Convert data from wide to long format for ANOVA
df_long <- melt(df, id.vars = "Mesh_Type", variable.name = "Control_Points", value.name = "Value")

# Subset the data for Poisson meshes and Aligned-Only
poisson_data <- subset(df_long, Mesh_Type == "Poisson")$Value
aligned_data <- subset(df_long, Mesh_Type == "Aligned-Only")$Value

# Step 1: Test for normality
# Shapiro-Wilk test for normality
shapiro.test(df_long$Value) #Data is normal 
leveneTest(df_long$Value, df_long$Mesh_Type)

# Step 2: Conduct the t-test
# Welch Two Sample t-test
t_test_result <- t.test(poisson_data, aligned_data)
print(t_test_result)

# Step 3: Standard Deviation 
sd_poisson <- sd(poisson_data)
print(sd_poisson)

sd_aligned <- sd(aligned_data)
print(sd_aligned)
