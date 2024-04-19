# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Compare eigenvalue scores 

# Author: James M. Mulqueeney 

# Date Last Modified: 06/02/2024 

# Poison Mesh Comparison with Manual Landmarking Data 

# Read in Libraries 
library(rgl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(geomorph)

#########################################################################################

# Read shape data & species data 
shape.data <- read.csv("E:/CTData/James Mulqueeney/Mammalian Data/Full Results/Manual Landmarking/Data/shape.data.322.csv")
species.data <- read.csv("E:/CTData/James Mulqueeney/Mammalian Data/Full Results/Manual Landmarking/Data/full_species_data.csv")

# Perform PCA on shape data 
PCA1 <- gm.prcomp(shape.data[,-1])

# Extract summary and write to a .txt file 
PCAsummary <- summary(PCA1)
write.table(PCAsummary$PC.summary, file = "pcasummary.txt")

# Extract eigenvalues from PCA results
eigenvalues <- as.numeric(t(PCAsummary$PC.summary[2,]) * 100)

#########################################################################################

# Aligned meshes 

# Read Kernel 40.0 eigenvalues data 
Kernel_40.0_eig_a <- read.csv("E:/CTData/James Mulqueeney/Papers/Write Up Papers/Paper 2- Comparison of Methods/Data/Results/Landmark Aligned/Kernel 40.0/eigenvalues.csv")

# Initialize empty vector for eigenvalues2
eigenvalues2 <- as.numeric(nrow(Kernel_40.0_eig_a))

# Calculate the eigenvalues iteratively
for (i in 1:nrow(Kernel_40.0_eig_a)) {
  if (i == 1) {
    eigenvalues2[i] <- Kernel_40.0_eig_a$cum..variability..in...[i]
  } else {
    eigenvalues2[i] <- Kernel_40.0_eig_a$cum..variability..in...[i] - Kernel_40.0_eig_a$cum..variability..in...[i - 1]
  }
}

# Print the resulting dataframe
print(eigenvalues2)

#########################################################################################

# Read Kernel 20.0 eigenvalues data 
Kernel_20.0_eig_a <- read.csv("E:/CTData/James Mulqueeney/Papers/Write Up Papers/Paper 2- Comparison of Methods/Data/Results/Landmark Aligned/Kernel 20.0/eigenvalues.csv")

# Initialize empty vector for eigenvalues3
eigenvalues3 <- as.numeric(nrow(Kernel_20.0_eig_a))

# Calculate the eigenvalues iteratively
for (i in 1:nrow(Kernel_20.0_eig_a)) {
  if (i == 1) {
    eigenvalues3[i] <- Kernel_20.0_eig_a$cum..variability..in...[i]
  } else {
    eigenvalues3[i] <- Kernel_20.0_eig_a$cum..variability..in...[i] - Kernel_20.0_eig_a$cum..variability..in...[i - 1]
  }
}

# Print the resulting dataframe
print(eigenvalues3)

#########################################################################################

# Read Kernel 10.0 eigenvalues data 
Kernel_10.0_eig_a <- read.csv("E:/CTData/James Mulqueeney/Papers/Write Up Papers/Paper 2- Comparison of Methods/Data/Results/Landmark Aligned/Kernel 10.0/eigenvalues.csv")

# Initialize empty vector for eigenvalues4
eigenvalues4 <- as.numeric(nrow(Kernel_10.0_eig_a))

# Calculate the eigenvalues iteratively
for (i in 1:nrow(Kernel_10.0_eig_a)) {
  if (i == 1) {
    eigenvalues4[i] <- Kernel_10.0_eig_a$cum..variability..in...[i]
  } else {
    eigenvalues4[i] <- Kernel_10.0_eig_a$cum..variability..in...[i] - Kernel_10.0_eig_a$cum..variability..in...[i - 1]
  }
}

# Print the resulting dataframe
print(eigenvalues4)

#########################################################################################

# Poisson meshes 

# Read Kernel 40.0 eigenvalues data 
Kernel_40.0_eig <- read.csv("E:/CTData/James Mulqueeney/Papers/Write Up Papers/Paper 2- Comparison of Methods/Data/Results/Poisson Meshes/Kernel 40.0/eigenvalues.csv")

# Initialize empty vector for eigenvalues2
eigenvalues5 <- as.numeric(nrow(Kernel_40.0_eig))

# Calculate the eigenvalues iteratively
for (i in 1:nrow(Kernel_40.0_eig)) {
  if (i == 1) {
    eigenvalues5[i] <- Kernel_40.0_eig$cum..variability..in...[i]
  } else {
    eigenvalues5[i] <- Kernel_40.0_eig$cum..variability..in...[i] - Kernel_40.0_eig$cum..variability..in...[i - 1]
  }
}

# Print the resulting dataframe
print(eigenvalues5)

#########################################################################################

# Read Kernel 20.0 eigenvalues data 
Kernel_20.0_eig <- read.csv("E:/CTData/James Mulqueeney/Papers/Write Up Papers/Paper 2- Comparison of Methods/Data/Results/Poisson Meshes/Kernel 20.0/eigenvalues.csv")

# Initialize empty vector for eigenvalues3
eigenvalues6 <- as.numeric(nrow(Kernel_20.0_eig))

# Calculate the eigenvalues iteratively
for (i in 1:nrow(Kernel_20.0_eig)) {
  if (i == 1) {
    eigenvalues6[i] <- Kernel_20.0_eig$cum..variability..in...[i]
  } else {
    eigenvalues6[i] <- Kernel_20.0_eig$cum..variability..in...[i] - Kernel_20.0_eig$cum..variability..in...[i - 1]
  }
}

# Print the resulting dataframe
print(eigenvalues6)

#########################################################################################

# Read Kernel 10.0 eigenvalues data 
Kernel_10.0_eig <- read.csv("E:/CTData/James Mulqueeney/Papers/Write Up Papers/Paper 2- Comparison of Methods/Data/Results/Poisson Meshes/Kernel 10.0/eigenvalues.csv")

# Initialize empty vector for eigenvalues4
eigenvalues7 <- as.numeric(nrow(Kernel_10.0_eig))

# Calculate the eigenvalues iteratively
for (i in 1:nrow(Kernel_10.0_eig)) {
  if (i == 1) {
    eigenvalues7[i] <- Kernel_10.0_eig$cum..variability..in...[i]
  } else {
    eigenvalues7[i] <- Kernel_10.0_eig$cum..variability..in...[i] - Kernel_10.0_eig$cum..variability..in...[i - 1]
  }
}

# Print the resulting dataframe
print(eigenvalues7)

########################################################################################

# Create a dataframe for eigenvalues
eigen_df <- data.frame(eigenvalues, eigenvalues2, eigenvalues3, eigenvalues4, eigenvalues5, eigenvalues6, eigenvalues7)

# Print or view the dataframe
print(eigen_df)

# Reshape the dataframe to long format
eigen_df_long <- reshape2::melt(eigen_df, id.vars = "eigenvalues")

# Plot
ggplot(eigen_df_long, aes(x = eigenvalues, y = value, color = variable)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, size=1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Manual Landmarking", y = "Deterministic Atlas Analysis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Plot
ggplot(eigen_df_long, aes(x = log(eigenvalues), y = log(value), color = variable)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(method = "lm", se = TRUE, size=1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Manual Landmarking", y = "Deterministic Atlas Analysis") +
  xlim(-17, 5)+ ylim(-17, 5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Calculate correlations
correlation_eigenvalues2 <- cor(eigenvalues, eigenvalues2)
correlation_eigenvalues3 <- cor(eigenvalues, eigenvalues3)
correlation_eigenvalues4 <- cor(eigenvalues, eigenvalues4)
correlation_eigenvalues5 <- cor(eigenvalues, eigenvalues5)
correlation_eigenvalues6 <- cor(eigenvalues, eigenvalues6)
correlation_eigenvalues7 <- cor(eigenvalues, eigenvalues7)

# Print correlations
print(correlation_eigenvalues2) # Aligned Only - k = 40
print(correlation_eigenvalues3) # Aligned Only - k = 20
print(correlation_eigenvalues4) # Aligned Only - k = 10
print(correlation_eigenvalues5) # Poisson = k = 40
print(correlation_eigenvalues6) # Poisson = k = 20
print(correlation_eigenvalues7) # Poisson = k = 10

########################################################################################

# Subset the eigenvalues variables to include only the first four PC axes 
eigen_df_4 <- data.frame(eigenvalues[1:4], eigenvalues2[1:4], eigenvalues3[1:4], eigenvalues4[1:4])

# Print or view the dataframe
print(eigen_df_4)

# Reshape the dataframe to long format
eigen_df_long_4 <- reshape2::melt(eigen_df_4, id.vars = "eigenvalues.1.4.")

# Plot
ggplot(eigen_df_long_4, aes(x = eigenvalues.1.4., y = value, color = variable)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, size=1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Manual Landmarking", y = "Deterministic Atlas Analysis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Calculate correlations
correlation_eigenvalues2 <- cor(eigenvalues[1:2], eigenvalues2[1:4])
correlation_eigenvalues3 <- cor(eigenvalues[1:2], eigenvalues3[1:4])
correlation_eigenvalues4 <- cor(eigenvalues[1:2], eigenvalues4[1:4])

# Print correlations
print(correlation_eigenvalues2)
print(correlation_eigenvalues3)
print(correlation_eigenvalues4)





