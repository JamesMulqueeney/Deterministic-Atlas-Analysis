# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Poisson Mesh - Interactive Plots

################################################################################

# Load required libraries
library(plotly)

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
  p1 <- plot_ly(data, x = ~PC1, y = ~PC2) %>%
    add_markers(text = ~Taxon, color = ~Order, colors = colortable$Order_colour, 
                symbol = ~Mesh.Type, symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
                marker = list(line = list(color = ifelse(data$Status == "Extinct", "gray", "black"), width = 0.5)))
  
  p2 <- plot_ly(data, x = ~PC3, y = ~PC4) %>%
    add_markers(text = ~Taxon, color = ~Order, colors = colortable$Order_colour, 
                symbol = ~Mesh.Type, symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
                marker = list(line = list(color = ifelse(data$Status == "Extinct", "gray", "black"), width = 0.5)))
  
  p3 <- plot_ly(data, x = ~PC1, y = ~PC2, z = ~PC3) %>%
    add_markers(text = ~Taxon, color = ~Order, colors = colortable$Order_colour, 
                symbol = ~Mesh.Type, symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
                marker = list(line = list(color = ifelse(data$Status == "Extinct", "gray", "black"), width = 0.5)))
  
  list(p1, p2, p3)
}

################################################################################

# Apply to the Poisson DAA results 

################################################################################

# Read in Poisson DAA results 
Kernel_40.0 <- read.csv(file.path(input_dir, "Data_A7-Poisson_k40_kpca.csv"))
Kernel_20.0 <- read.csv(file.path(input_dir, "Data_A8-Poisson_k20_kpca.csv"))
Kernel_10.0 <- read.csv(file.path(input_dir, "Data_A9-Poisson_k10_kpca.csv"))

# Generate the plots 
plots_40 <- generate_plots(Kernel_40.0, "Kernel 40.0")
plots_20 <- generate_plots(Kernel_20.0, "Kernel 20.0")
plots_10 <- generate_plots(Kernel_10.0, "Kernel 10.0")

# Display plots
plots_40
plots_20
plots_10

