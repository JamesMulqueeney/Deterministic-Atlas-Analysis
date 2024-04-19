# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies on the mammalian crania 

# Produce Principal Component Plots 

# Author: James M. Mulqueeney 

# Date Last Modified: 13/03/2024

# Manual Landmarking Data Only 

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
shape.data <- read.csv("path/to/input/Data_S3-Shape_Data_322.csv")

# Species data (Details of Taxonomy etc.)
species.data <- read.csv("path/to/input/Data_S1-Specimen_Details.csv")

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
colortable <- read.csv("path/to/input/Data_S17-New_Order_Colors.csv")
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

# Initial Plot of Manual Landmarking Data with the Figure Legend 

# Plot the Data with Legend 
g1 <- ggplot(pcscores, aes(x=Comp1, y=Comp2, shape=Mesh.Type,fill=Order,color=Status)) +
  geom_point(size=1.5,stroke=.4)+
  scale_shape_manual(name = "", values=c(21, 22,23,25, 24))+
  scale_color_manual(name = "", values = c("black","grey60"))+
  scale_fill_manual(values = neworder)+
  labs(
    x = paste0("Principal Component 1 (", signif((PCAsummary$PC.summary[2,1]*100),3), "%)"),
    y = paste0("Principal Component 2 (", signif((PCAsummary$PC.summary[2,2]*100),3),"%)")
  )+
  scale_x_continuous(breaks = seq(-0.1,0.5,by=.1))+
  scale_y_continuous(breaks = seq(-0.2,0.3,by=.1))+
  theme_bw()+
  coord_fixed(ratio=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=9),
        legend.text=element_text(size=13));g1
g1<-g1+guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(fill = guide_legend(override.aes=list(shape=shapevector, color="white", size=5)))+
  guides(shape = guide_legend(override.aes=list(size=5)))
g1

# Save the Plot g1

# Specify the file path including the desired directory and filename
g1_file_path <- "path\\to\\output\\filename.png"

# Use ggsave to save the plot to the specified file path
ggsave(g1_file_path, plot = g1, device = "png", width = 10, height = 6, units = "in", dpi = 600)

#########################################################################################

##############################################
# Manual Landmarking plots without the legend
##############################################

############
# PC1 & PC2
############

g2 <- ggplot(pcscores, aes(x = Comp1, y = Comp2, shape = Mesh.Type, fill = Order, color = Status)) +
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

g2

# Save the Plot g2

# Specify the file path including the desired directory and filename
g2_file_path <- "path\\to\\output\\filename.png"

# Use ggsave to save the plot to the specified file path
ggsave(g2_file_path, plot = g2, device = "png", width = 10, height = 6, units = "in", dpi = 600)


#########################################################################################

############
# PC3 & PC4
############

g3 <- ggplot(pcscores, aes(x = Comp3, y = Comp4, shape = Mesh.Type, fill = Order, color = Status)) +
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

g3

# Save the Plot g3

# Specify the file path including the desired directory and filename
g3_file_path <- "path\\to\\output\\filename.png"

# Use ggsave to save the plot to the specified file path
ggsave(g3_file_path, plot = g3, device = "png", width = 10, height = 6, units = "in", dpi = 600)

