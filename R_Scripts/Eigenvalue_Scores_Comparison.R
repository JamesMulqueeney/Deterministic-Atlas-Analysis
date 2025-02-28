# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Eigenvalue Comparisons 

################################################################################

# Load required libraries
library(geomorph)
library(ggplot2)

################################################################################

# Define input directory (change accordingly)
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data"

################################################################################

# Process the manual landmarking data 

################################################################################

# Read species data
species.data <- read.csv(file.path(input_dir, "Data_A20-Specimen_Details.csv"))

# Read in manual landmarking data 
shape.data <- read.csv(file.path(input_dir, "Data_A21-Shape_Data_322.csv"))

# Perform PCA on shape data 
PCA1 <- gm.prcomp(shape.data[,-1])

# Extract summary and write to a .txt file 
PCAsummary <- summary(PCA1)
write.table(PCAsummary$PC.summary, file = "pcasummary.txt")

# Extract eigenvalues from PCA results
eigenvalues <- as.numeric(t(PCAsummary$PC.summary[2,]) * 100)

################################################################################

# Aligned-only Mesh data comparison 

################################################################################

# Read Kernel 40.0 eigenvalues data 
Kernel_40.0_eig_a <- read.csv(file.path(input_dir, "Data_A12-Aligned_Only_k40_eigenvalues.csv"))

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

################################################################################

# Read Kernel 20.0 eigenvalues data 
Kernel_20.0_eig_a <- read.csv(file.path(input_dir, "Data_A13-Aligned_Only_k20_eigenvalues.csv"))

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

################################################################################

# Read Kernel 10.0 eigenvalues data 
Kernel_10.0_eig_a <- read.csv(file.path(input_dir, "Data_A14-Aligned_Only_k10_eigenvalues.csv"))

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

################################################################################

# Poisson meshes 

################################################################################

# Read Kernel 40.0 eigenvalues data 
Kernel_40.0_eig <- read.csv(file.path(input_dir, "Data_A15-Poisson_k40_eigenvalues.csv"))

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

################################################################################

# Read Kernel 20.0 eigenvalues data 
Kernel_20.0_eig <- read.csv(file.path(input_dir, "Data_A16-Poisson_k20_eigenvalues.csv"))

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

################################################################################

# Read Kernel 40.0 eigenvalues data 
Kernel_10.0_eig <- read.csv(file.path(input_dir, "Data_A17-Poisson_k10_eigenvalues.csv"))

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

################################################################################

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
correlation_eigenvalues2 <- lm(eigenvalues ~ eigenvalues2)
correlation_eigenvalues3 <- lm(eigenvalues ~ eigenvalues3)
correlation_eigenvalues4 <- lm(eigenvalues ~ eigenvalues4)
correlation_eigenvalues5 <- lm(eigenvalues ~ eigenvalues5)
correlation_eigenvalues6 <- lm(eigenvalues ~ eigenvalues6)
correlation_eigenvalues7 <- lm(eigenvalues ~ eigenvalues7)

# Print correlations
summary(correlation_eigenvalues2) # Aligned Only - k = 40 # 0.9421
summary(correlation_eigenvalues3) # Aligned Only - k = 20 # 0.9617
summary(correlation_eigenvalues4) # Aligned Only - k = 10 # 0.9753
summary(correlation_eigenvalues5) # Poisson = k = 40 # 0.9934
summary(correlation_eigenvalues6) # Poisson = k = 20 # 0.9775
summary(correlation_eigenvalues7) # Poisson = k = 10 # 0.9619

########################################################################################

# Subset the eigenvalues variables to include only the first four PC axes 
eigen_df_4 <- data.frame(eigenvalues[1:10], eigenvalues2[1:10], eigenvalues3[1:10], eigenvalues4[1:10])

# Print or view the dataframe
print(eigen_df_4)

# Reshape the dataframe to long format
eigen_df_long_4 <- reshape2::melt(eigen_df_4, id.vars = "eigenvalues.1.10.")

# Plot
ggplot(eigen_df_long_4, aes(x = eigenvalues.1.10., y = value, color = variable)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, size=1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  labs(x = "Manual Landmarking", y = "Deterministic Atlas Analysis") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Calculate correlations
correlation_eigenvalues2.2 <- lm(eigenvalues[1:10] ~ eigenvalues2[1:10])
correlation_eigenvalues3.2 <- lm(eigenvalues[1:10] ~ eigenvalues3[1:10])
correlation_eigenvalues4.2 <- lm(eigenvalues[1:10] ~ eigenvalues4[1:10])
correlation_eigenvalues5.2 <- lm(eigenvalues[1:10] ~ eigenvalues5[1:10])
correlation_eigenvalues6.2 <- lm(eigenvalues[1:10] ~ eigenvalues6[1:10])
correlation_eigenvalues7.2 <- lm(eigenvalues[1:10] ~ eigenvalues7[1:10])

# Print correlations
summary(correlation_eigenvalues2.2) # Aligned Only - k = 40 # 0.9587
summary(correlation_eigenvalues3.2) # Aligned Only - k = 20 # 0.9546
summary(correlation_eigenvalues4.2) # Aligned Only - k = 10 # 0.9665
summary(correlation_eigenvalues5.2) # Poisson = k = 40 # 0.9931
summary(correlation_eigenvalues6.2) # Poisson = k = 20 # 0.9917
summary(correlation_eigenvalues7.2) # Poisson = k = 10 # 0.9921

