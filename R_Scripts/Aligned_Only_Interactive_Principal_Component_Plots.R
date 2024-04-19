# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Produce Principal Component Plots 

# Author: James M. Mulqueeney 

# Date Last Modified: 13/03/2024

# Aligned Only Comparison with Manual Landmarking Data 

# Read in Libraries 
library(plotly)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)

#########################################################################################

# Read in colour/ shape of plot data
colortable <- read.csv("path/to/input/new_order_colors2.csv")
names(colortable$Order_colour)<-colortable$Order

# Re-order the Orders (Needs working on)
phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea","Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", "Astrapotheria","Perissodactyla", "Creodonta", "Carnivora", "Pholidota")

#######################################
# Kernel 40.0 (Control points = 45 )
#######################################

# Read in the data 
Kernel_40.0 <- read.csv("path/to/input/Data_S5-Aligned_Only_k40_kpca.csv")

# Plot in 2D 

# PC1 & PC2
plot1 <- plot_ly(Kernel_40.0, x = ~PC1, y = ~PC2) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_40.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot1

# PC3 & PC4
plot2 <- plot_ly(Kernel_40.0, x = ~PC3, y = ~PC4) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_40.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot2

# Plot in 3D 

# PC1 - PC3
plot3 <- plot_ly(Kernel_40.0, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Typer,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_40.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot3

#########################################################################################

#######################################
# Kernel 20.0 (Control points = 270 )
#######################################

# Read in the data 
Kernel_20.0 <- read.csv("path/to/input/Data_S6-Aligned_Only_k20_kpca.csv")

# Plot in 2D 

# PC1 & PC2
plot4 <- plot_ly(Kernel_20.0, x = ~PC1, y = ~PC2) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_20.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot4

# PC3 & PC4
plot5 <- plot_ly(Kernel_20.0, x = ~PC3, y = ~PC4) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_20.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot5

# Plot in 3D 

# PC1 - PC3
plot6 <- plot_ly(Kernel_20.0, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_20.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot6

#########################################################################################

#######################################
# Kernel 10.0 (Control points = 1782 )
#######################################

# Read in the data 
Kernel_10.0 <- read.csv("path/to/input/Data_S7-Aligned_Only_k10_kpca.csv")

# Plot in 2D 

# PC1 & PC2
plot7 <- plot_ly(Kernel_10.0, x = ~PC1, y = ~PC2) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_10.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot7

# PC3 & PC4
plot8 <- plot_ly(Kernel_10.0, x = ~PC3, y = ~PC4) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_10.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot8

# Plot in 3D 

# PC1 - PC3
plot9 <- plot_ly(Kernel_10.0, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  add_markers(text = ~Taxon, color = ~Order,colors = colortable$Order_colour,
              symbol = ~Mesh.Type,symbols = c("circle", "square", "diamond", "triangle-down", "triangle-up"),
              marker = list(line = list(color = ifelse(Kernel_10.0$Status == "Extinct", "gray", "black"), width = 0.5)),
              legendgroup = ~Order) %>%
  layout(legend = list(traceorder = "normal",tracegroupgap = 10,itemsizing = "constant",itemorder = "phyloseq",itemarray = "phyloseq"))

plot9
