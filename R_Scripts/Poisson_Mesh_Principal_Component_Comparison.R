# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Produce Principal Component Plots 

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

################################################
# Manual Landmarking (from Goswami et al. 2022)
################################################

# Read in shape data & species data 

# Shape data (Landmark Data)
shape.data <- read.csv("path/to/input/shape.data.322.csv")

# Species data (Details of Taxonomy etc.)
species.data <- read.csv("path/to/input/full_species_data.csv")

#########################################################################################

# Perform PCA on shape data 
PCA1 <- gm.prcomp(shape.data[,-1])

# Extract summary and write to a .txt file 
PCAsummary<-summary(PCA1)
write.table(PCAsummary$PC.summary, file = "pcasummary.txt")

# Extract the PCA values 
pcscores<-PCA1$x

# Join PCA and Species data together 
pcscores<-cbind(species.data,pcscores) 

# Assign colours 
colortable <- read.csv("E:/CTData/James Mulqueeney/Mammalian Data/Full Results/Coding Files/new_order_colors2.csv")
names(colortable$Order_colour)<-colortable$Order
#names(colortable$order_color3) <- colortable$Order

# Assign to the colours to the correct Orders 
phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea","Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", "Astrapotheria","Perissodactyla", "Creodonta", "Carnivora", "Pholidota")
neworder <- colortable$Order_colour
# neworder <- colortable$order_color3
neworder <- colortable$Order_colour[order(factor(colortable$Order_colour, levels = phyloseq))]
#pcscores$Order <- factor(pcscores$Original_Order, levels = phyloseq)
#pcscores <- pcscores %>% mutate(order2 = Order)

# Assign the shape to the correct Orders  
names(colortable$ordershape) <- colortable$Order
shapevector <- colortable$ordershape
shapevector <- colortable$ordershape[order(factor(colortable$ordershape, levels = phyloseq))]
shapevector <- as.vector(shapevector)
# pcscores$order2[which(pcscores$Extant_Extinct=="Extinct")]<-NA

#########################################################################################

############
# PC1 & PC2
############

g1 <- ggplot(pcscores, aes(x = Comp1, y = Comp2, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 1 (", signif((PCAsummary$PC.summary[2,1]*100),3), "%)"),
    y = paste0("Principal Component 2 (", signif((PCAsummary$PC.summary[2,2]*100),3),"%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g1

# Save the Plot g1

# Specify the file path including the desired directory and filename
g1_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Manual Landmarking\\Manual_Landmarking_PC1-2.png"

# Use ggsave to save the plot to the specified file path
ggsave(g1_file_path, plot = g1, device = "png", width = 10, height = 6, units = "in", dpi = 600)

############
# PC3 & PC4
############

g2 <- ggplot(pcscores, aes(x = Comp3, y = Comp4, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 3 (", signif((PCAsummary$PC.summary[2,3]*100),3), "%)"),
    y = paste0("Principal Component 4 (", signif((PCAsummary$PC.summary[2,4]*100),3),"%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g2

# Save the Plot g2

# Specify the file path including the desired directory and filename
g2_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Manual Landmarking\\Manual_Landmarking_PC3-4.png"

# Use ggsave to save the plot to the specified file path
ggsave(g2_file_path, plot = g2, device = "png", width = 10, height = 6, units = "in", dpi = 600)

#########################################################################################

#####################################################
# Deterministic Atlas Analysis (DAA) : Aligned-Only 
#####################################################

#########################################################################################

#######################################
# Kernel 40.0 (Control points = 45 )
#######################################

# Read in the data 
Kernel_40.0 <- read.csv("path/to/input//Poisson Meshes/Kernel 40.0/kpca.csv")

############
# PC1 & PC2
############

g3 <- ggplot(Kernel_40.0, aes(x=PC1, y=PC2, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 1 (39.17%)"),
    y = paste0("Principal Component 2 (18.95%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g3

# Save the Plot g3

# Specify the file path including the desired directory and filename
g3_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Deterministic Atlas Analysis\\Kernel 40.0\\Poisson Meshes\\Kernel_40.0_PC1-2.png"

# Use ggsave to save the plot to the specified file path
ggsave(g3_file_path, plot = g3, device = "png", width = 10, height = 6, units = "in", dpi = 600)

############
# PC3 & PC4
############

g4 <- ggplot(Kernel_40.0, aes(x=PC3, y=PC4, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 3 (12.58%)"),
    y = paste0("Principal Component 4 (6.42%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g4

# Save the Plot g4

# Specify the file path including the desired directory and filename
g4_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Deterministic Atlas Analysis\\Kernel 40.0\\Poisson Meshes\\Kernel_40.0_PC3-4.png"

# Use ggsave to save the plot to the specified file path
ggsave(g4_file_path, plot = g4, device = "png", width = 10, height = 6, units = "in", dpi = 600)

#########################################################################################

#######################################
# Kernel 20.0 (Control points = 270 )
#######################################

# Read in the data 
Kernel_20.0 <- read.csv("path/to/input/Poisson Meshes/Kernel 20.0/kpca.csv")

############
# PC1 & PC2
############

g5 <- ggplot(Kernel_20.0, aes(x=PC1, y=PC2, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 1 (25.62%)"),
    y = paste0("Principal Component 2 (13.76%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g5

# Save the Plot g5

# Specify the file path including the desired directory and filename
g5_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Deterministic Atlas Analysis\\Kernel 20.0\\Poisson Meshes\\Kernel_20.0_PC1-2.png"

# Use ggsave to save the plot to the specified file path
ggsave(g5_file_path, plot = g5, device = "png", width = 10, height = 6, units = "in", dpi = 600)

############
# PC3 & PC4
############

g6 <- ggplot(Kernel_20.0, aes(x=PC3, y=PC4, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 3 (8.88%)"),
    y = paste0("Principal Component 4 (5.83%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g6

# Save the Plot g6

# Specify the file path including the desired directory and filename
g6_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Deterministic Atlas Analysis\\Kernel 20.0\\Poisson Meshes\\Kernel_20.0_PC3-4.png"

# Use ggsave to save the plot to the specified file path
ggsave(g6_file_path, plot = g6, device = "png", width = 10, height = 6, units = "in", dpi = 600)

#########################################################################################

#######################################
# Kernel 10.0 (Control points = 1782 )
#######################################

# Read in the data 
Kernel_10.0 <- read.csv("path/to/input/Poisson Meshes/Kernel 10.0/kpca.csv")

############
# PC1 & PC2
############

g7 <- ggplot(Kernel_10.0, aes(x=PC1, y=PC2, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 1 (19.01%)"),
    y = paste0("Principal Component 2 (9.65%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g7

# Save the Plot g7

# Specify the file path including the desired directory and filename
g7_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Deterministic Atlas Analysis\\Kernel 10.0\\Poisson Meshes\\Kernel_10.0_PC1-2.png"

# Use ggsave to save the plot to the specified file path
ggsave(g7_file_path, plot = g7, device = "png", width = 10, height = 6, units = "in", dpi = 600)

############
# PC3 & PC4
############

g8 <- ggplot(Kernel_10.0, aes(x=PC3, y=PC4, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = paste0("Principal Component 3 (6.75%)"),
    y = paste0("Principal Component 4 (5.49%)")
  )+
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 13)
  )

g8

# Save the Plot g8

# Specify the file path including the desired directory and filename
g8_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Deterministic Atlas Analysis\\Kernel 10.0\\Poisson Meshes\\Kernel_10.0_PC3-4.png"

# Use ggsave to save the plot to the specified file path
ggsave(g8_file_path, plot = g8, device = "png", width = 10, height = 6, units = "in", dpi = 600)

#########################################################################################

###################################
# Combined Figures (PC Comparisons)
###################################

###########
# PC1 & PC2 
###########

# Combined Plot 
g9 <- plot_grid(
  g1, g3, g5, g7,
  nrow = 2, byrow = TRUE,
  labels = c("A", "B", "C", "D"),
  align="hv")

g9

# Save the Plot g9

# Specify the file path including the desired directory and filename
g9_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Combined PC Plots\\Poisson Meshes\\Combined_PC1_PC2.png"

# Use ggsave to save the plot to the specified file path
ggsave(g9_file_path, plot = g9, device = "png", width = 10, height = 6, units = "in", dpi = 600)

###########
# PC3 & PC4 
###########

# Combined Plot 
g10 <- plot_grid(
  g2, g4, g6, g8,
  nrow = 2, byrow = TRUE,
  labels = c("A", "B", "C", "D"),
  align="hv")

g10

# Save the Plot g10

# Specify the file path including the desired directory and filename
g10_file_path <- "E:\\CTData\\James Mulqueeney\\Papers\\Write Up Papers\\Paper 2- Comparison of Methods\\Figures\\Combined PC Plots\\Poisson Meshes\\Combined_PC3_PC4.png"

# Use ggsave to save the plot to the specified file path
ggsave(g10_file_path, plot = g10, device = "png", width = 10, height = 6, units = "in", dpi = 600)

#########################################################################################
