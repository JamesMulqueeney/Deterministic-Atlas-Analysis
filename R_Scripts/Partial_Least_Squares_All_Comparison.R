# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Partial Least Sqaures Comparison 

################################################################################

# Load required libraries 
library(Morpho)
library(abind)
library(Rvcg)
library(dplyr)
library(geomorph)
library(ggplot2)
library(purrr)

################################################################################

# Define input directory
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/"

################################################################################

# Provide plotting information 

################################################################################

# Read species data
species.data <- read.csv(file.path(input_dir, "Data_A20-Specimen_Details.csv"))
colortable <- read.csv(file.path(input_dir, "Data_A23-New_Order_Colors.csv"))
names(colortable$Order_colour) <- colortable$Order

# Define order levels
phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", 
              "Macroscelidea", "Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", 
              "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", 
              "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", 
              "Astrapotheria", "Perissodactyla", "Creodonta", "Carnivora", "Pholidota")

# Assign colors and shapes
neworder <- colortable$Order_colour[order(factor(colortable$Order_colour, levels = phyloseq))]
names(colortable$ordershape) <- colortable$Order
shapevector <- as.vector(colortable$ordershape[order(factor(colortable$ordershape, levels = phyloseq))])

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

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result1 <- two.b.pls(merged_data1[, comps], merged_data1[, pcs])
summary(pls_result1)

# Visualize the results
PLS1<- plot(pls_result1, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

# Make a Better Plot

# Extract scores from the PLS analysis
scores1.1 <- pls_result1$XScores[, 1]
scores2.1 <- pls_result1$YScores[, 1]

# Create a dataframe with scores for plotting
pls_df1 <- data.frame(scores1.1, scores2.1)

# Add additional variables to the dataframe if needed
# For example, assuming Kernel_40.0 contains additional variables:
pls_df1$Mesh.Type <- merged_data1$Mesh.Type.x
pls_df1$Order <- merged_data1$Order.x
pls_df1$Status <- merged_data1$Status.x

# Create the ggplot plot
PLS1 <- ggplot(pls_df1, aes(x = scores1.1, y = scores2.1, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size=1) +  # Add a single regression line for the entire dataset
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  xlim(-0.5, 0.2) + ylim(-0.2, 0.15)+
  labs(
    x = paste0("Manual Landmarking"),
    y = paste0("Deterministic Atlas Analysis")
  ) +
  theme_bw() +
  coord_fixed(ratio = 2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 13)
  )

# Print the plot
print(PLS1)

################################################################################

# Aligned_Only: Manual vs Kernel 20.0 (Control points = 270 )

################################################################################

# Merge the data into a new dataframe 
merged_data2 <- merge(pcscores, Kernel_20.0, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result2 <- two.b.pls(merged_data2[, comps], merged_data2[, pcs])
summary(pls_result2)

# Visualize the results
PLS2<- plot(pls_result2, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

# Make a Better Plot

# Extract scores from the PLS analysis
scores1.2 <- pls_result2$XScores[, 1]
scores2.2 <- pls_result2$YScores[, 1]

# Create a dataframe with scores for plotting
pls_df2 <- data.frame(scores1.2, scores2.2)

# Add additional variables to the dataframe if needed
# For example, assuming Kernel_40.0 contains additional variables:
pls_df2$Mesh.Type <- merged_data1$Mesh.Type.x
pls_df2$Order <- merged_data1$Order.x
pls_df2$Status <- merged_data1$Status.x

# Create the ggplot plot
PLS2 <- ggplot(pls_df2, aes(x = scores1.2, y = scores2.2, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size=1) + 
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  xlim(-0.5, 0.2) + ylim(-0.2, 0.15)+
  labs(
    x = paste0("Manual Landmarking"),
    y = paste0("Deterministic Atlas Analysis")
  ) +
  theme_bw() +
  coord_fixed(ratio = 2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 13)
  )

# Print the plot
print(PLS2)

################################################################################

# Aligned_Only: Manual vs Kernel 10.0 (Control points = 1782 )

################################################################################

# Merge the data into a new dataframe 
merged_data3 <- merge(pcscores, Kernel_10.0, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result3 <- two.b.pls(merged_data3[, comps], merged_data3[, pcs])
summary(pls_result3)

# Visualize the results
PLS3<- plot(pls_result3, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

# Make a Better Plot

# Extract scores from the PLS analysis
scores1.3 <- pls_result3$XScores[, 1]
scores2.3 <- pls_result3$YScores[, 1]

# Create a dataframe with scores for plotting
pls_df3 <- data.frame(scores1.3, scores2.3)

# Add additional variables to the dataframe if needed
# For example, assuming Kernel_40.0 contains additional variables:
pls_df3$Mesh.Type <- merged_data1$Mesh.Type.x
pls_df3$Order <- merged_data1$Order.x
pls_df3$Status <- merged_data1$Status.x

# Create the ggplot plot
PLS3 <- ggplot(pls_df3, aes(x = scores1.3, y = scores2.3, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size=1) + 
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  xlim(-0.5, 0.2) + ylim(-0.2, 0.15)+
  labs(
    x = paste0("Manual Landmarking"),
    y = paste0("Deterministic Atlas Analysis")
  ) +
  theme_bw() +
  coord_fixed(ratio = 2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 13)
  )

# Print the plot
print(PLS3)

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

# Manual vs Kernel 40.0 

################################################################################

# Merge the data into a new dataframe 
merged_data4 <- merge(pcscores, Kernel_40.0_P, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result4 <- two.b.pls(merged_data4[, comps], merged_data4[, pcs])
summary(pls_result4)

# Visualize the results
PLS4<- plot(pls_result4, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

# Make a Better Plot

# Extract scores from the PLS analysis
scores1.4 <- pls_result4$XScores[, 1]
scores2.4 <- pls_result4$YScores[, 1]

# Create a dataframe with scores for plotting
pls_df4 <- data.frame(scores1.4, scores2.4)

# Add additional variables to the dataframe if needed
# For example, assuming Kernel_40.0 contains additional variables:
pls_df4$Mesh.Type <- merged_data1$Mesh.Type.x
pls_df4$Order <- merged_data1$Order.x
pls_df4$Status <- merged_data1$Status.x

# Create the ggplot plot
PLS4 <- ggplot(pls_df4, aes(x = scores1.4, y = scores2.4, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size=1) + 
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  xlim(-0.5, 0.2) + ylim(-0.2, 0.15)+
  labs(
    x = paste0("Manual Landmarking"),
    y = paste0("Deterministic Atlas Analysis")
  ) +
  theme_bw() +
  coord_fixed(ratio = 2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 13)
  )

# Print the plot
print(PLS4)

################################################################################

# Manual vs Kernel 20.0 

################################################################################

# Merge the data into a new dataframe 
merged_data5 <- merge(pcscores, Kernel_20.0_P, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result5 <- two.b.pls(merged_data5[, comps], merged_data5[, pcs])
summary(pls_result5)

# Visualize the results
PLS5<- plot(pls_result5, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

# Make a Better Plot

# Extract scores from the PLS analysis
scores1.5 <- pls_result5$XScores[, 1]
scores2.5 <- pls_result5$YScores[, 1]

# Create a dataframe with scores for plotting
pls_df5 <- data.frame(scores1.5, scores2.5)

# Add additional variables to the dataframe if needed
# For example, assuming Kernel_40.0 contains additional variables:
pls_df5$Mesh.Type <- merged_data1$Mesh.Type.x
pls_df5$Order <- merged_data1$Order.x
pls_df5$Status <- merged_data1$Status.x

# Create the ggplot plot
PLS5 <- ggplot(pls_df5, aes(x = scores1.5, y = scores2.5, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size=1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  xlim(-0.5, 0.2) + ylim(-0.2, 0.15)+
  labs(
    x = paste0("Manual Landmarking"),
    y = paste0("Deterministic Atlas Analysis")
  ) +
  theme_bw() +
  coord_fixed(ratio = 2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 13)
  )

# Print the plot
print(PLS5)

###############################################################################

# Manual vs Kernel 10.0 

################################################################################

# Merge the data into a new dataframe 
merged_data6 <- merge(pcscores, Kernel_10.0_P, by = "Tip_Label")

# Perform partial least squares analysis
pcs <- paste0("PC", 1:321)
comps <- paste0("Comp", 1:321)

pls_result6 <- two.b.pls(merged_data6[, comps], merged_data6[, pcs])
summary(pls_result6)

# Visualize the results
PLS6 <- plot(pls_result6, type = "p", pch=19, cex=1, xlim=c(-0.5,0.2), ylim=c(-0.2,0.15))

# Make a Better Plot

# Extract scores from the PLS analysis
scores1.6 <- pls_result6$XScores[, 1]
scores2.6 <- pls_result6$YScores[, 1]

# Create a dataframe with scores for plotting
pls_df6 <- data.frame(scores1.6, scores2.6)

# Add additional variables to the dataframe if needed
# For example, assuming Kernel_40.0 contains additional variables:
pls_df6$Mesh.Type <- merged_data1$Mesh.Type.x
pls_df6$Order <- merged_data1$Order.x
pls_df6$Status <- merged_data1$Status.x

# Create the ggplot plot
PLS6 <- ggplot(pls_df6, aes(x = scores1.6, y = scores2.6, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size=1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  xlim(-0.5, 0.2) + ylim(-0.2, 0.15)+
  labs(
    x = paste0("Manual Landmarking"),
    y = paste0("Deterministic Atlas Analysis")
  ) +
  theme_bw() +
  coord_fixed(ratio = 2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 13)
  )

# Print the plot
print(PLS6)

#########################################################################################

# For ggplot combined

################################################################################

# Combined Plot 
p1 <- plot_grid(
  PLS1,PLS2,PLS3,PLS4,PLS5,PLS6,
  nrow = 2, byrow = TRUE,
  labels = c("A", "B", "C", "D", "E", "F"),
  align="hv")

p1

