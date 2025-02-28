# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Manual Landmarking Data 

################################################################################

# Load required libraries
library(plotly)
library(geomorph)

################################################################################

# Define input directory (change accordingly)
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data"

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

# Define a function to make PCA plots 

################################################################################

plot_pca <- function(data, xcomp, ycomp, var_x, var_y) {
  xlabel <- paste0("Principal Component ", gsub("PC", "", xcomp), " (", var_x, "%)")
  ylabel <- paste0("Principal Component ", gsub("PC", "", ycomp), " (", var_y, "%)")
  
  p <- ggplot(data, aes_string(x = xcomp, y = ycomp, shape = "Mesh.Type", fill = "Order", color = "Status")) +
    geom_point(size = 2.5, stroke = 0.4) +
    scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
    scale_color_manual(name = "", values = c("black", "grey60")) +
    scale_fill_manual(values = neworder) +
    labs(x = xlabel, y = ylabel) +
    theme_bw() +
    coord_fixed(ratio = 1) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      legend.text = element_text(size = 13)
    )
  
  return(p)  # Returns the plot instead of saving it
}

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

# Apply plot function 

# PC1 & PC2
p1 <- plot_pca(pcscores, "Comp1", "Comp2", 
               paste0(signif((PCAsummary$PC.summary[2,1]*100), 3)),  
               paste0(signif((PCAsummary$PC.summary[2,2]*100), 3)))

# PC3 & PC4
p2 <- plot_pca(pcscores, "Comp3", "Comp4", 
               paste0(signif((PCAsummary$PC.summary[2,3]*100), 3)),  
               paste0(signif((PCAsummary$PC.summary[2,4]*100), 3)))

################################################################################

# Combine plots 

################################################################################

# PC1 & PC2 Combined 
p3 <- plot_grid(
  p1, p2,
  nrow = 1, byrow = TRUE,
  labels = c("A", "B"),
  align="hv")

p3
