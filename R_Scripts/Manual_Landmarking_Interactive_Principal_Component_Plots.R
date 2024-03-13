# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Produce Principal Component Plots 

# Author: James M. Mulqueeney 

# Date Last Modified: 13/03/2024 

# Manual landmarking data only  

# Read in Libraries 
library(rgl)
library(plotly)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)

#########################################################################################

# Read in shape data & species data 

# Read in species data
Species.Data<-read.csv("path/to/input/full_species_data.csv")

# Manual landmarking results

# Read in shape data 
manual.shape.data<-read.csv("path/to/input/shape.data.322.csv")

# Perform PCA
PCA1<-prcomp(manual.shape.data[,-1])
summary(PCA1)

# Extract PC for each species
PCA_Values<-PCA1$x

# Combine PC values with the species data 
combined_data <- cbind(Species.Data, PCA_Values)

# Plot the data

# Read in colour/ shape of plot data
colortable <- read.csv("E:/CTData/James Mulqueeney/Mammalian Data/Full Results/Coding Files/new_order_colors2.csv")
names(colortable$Order_colour)<-colortable$Order

# Re-order the Orders (Needs working on)
phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea","Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", "Astrapotheria","Perissodactyla", "Creodonta", "Carnivora", "Pholidota")

#########################################################################################

# 2D Plot 

# PC1 & PC2 
plot <- plot_ly(combined_data, x = ~PC1, y = ~PC2) %>%
  add_markers(text = ~Tip_Label, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(combined_data$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot

# PC3 & PC4 
plot2 <- plot_ly(combined_data, x = ~PC3, y = ~PC4) %>%
  add_markers(text = ~Tip_Label, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(combined_data$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot2

# 3D Plot 

# PC1 - PC3
plot3 <- plot_ly(combined_data, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  add_markers(text = ~Tip_Label, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(combined_data$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot3



