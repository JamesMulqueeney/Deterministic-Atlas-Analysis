# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Manual Landmarking Data - Interactive Plots

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

# Define a function to make interactive PCA plots 

################################################################################

# Function to generate 2D and 3D plots
generate_plots <- function(data, kernel_name) {
  p1 <- plot_ly(data, x = ~Comp1, y = ~Comp2) %>%
    add_markers(text = ~Taxon, color = ~Order, colors = colortable$Order_colour, 
                symbol = ~Mesh.Type, symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
                marker = list(line = list(color = ifelse(data$Status == "Extinct", "gray", "black"), width = 0.5)))
  
  p2 <- plot_ly(data, x = ~Comp1, y = ~Comp2) %>%
    add_markers(text = ~Taxon, color = ~Order, colors = colortable$Order_colour, 
                symbol = ~Mesh.Type, symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
                marker = list(line = list(color = ifelse(data$Status == "Extinct", "gray", "black"), width = 0.5)))
  
  p3 <- plot_ly(data, x = ~Comp1, y = ~Comp2, z = ~Comp3) %>%
    add_markers(text = ~Taxon, color = ~Order, colors = colortable$Order_colour, 
                symbol = ~Mesh.Type, symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
                marker = list(line = list(color = ifelse(data$Status == "Extinct", "gray", "black"), width = 0.5)))
  
  list(p1, p2, p3)
}

################################################################################

# Apply to Manual landmarking results 

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

# Manual landmarking plot
pcscores_plot <- generate_plots(pcscores, "pcscores")
pcscores_plot
