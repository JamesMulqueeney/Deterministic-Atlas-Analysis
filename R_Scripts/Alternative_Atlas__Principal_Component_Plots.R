# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Alternative Atlas Comparison 

################################################################################

# Load required libraries
library(ggplot2)

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

# Apply to each Atlas result  

################################################################################

# Read in the data 
A_Binturong_Atlas <-  read.csv(file.path(input_dir, "Data_A1_A_binturong_Atlas_Results.csv"))
C_calvus_Atlas <- read.csv(file.path(input_dir, "Data_A2_C_calvus_Atlas_Results.csv"))
S_morckhoviensis_Atlas <- read.csv(file.path(input_dir, "Data_A3_S_morckhoviensis_Atals_Results.csv"))

# Apply plot function 
p1 <- plot_pca(A_Binturong_Atlas, "PC1", "PC2", 25.62, 13.76) # A_Binturong_Atlas - PC1-PC2
p2 <- plot_pca(A_Binturong_Atlas, "PC3", "PC4", 8.88, 5.83) # A_Binturong_Atlas- PC3-PC4
p3 <- plot_pca(C_calvus_Atlas, "PC1", "PC2", 22.41, 14.56) # C_calvus_Atlas - PC1-PC2
p4 <- plot_pca(C_calvus_Atlas, "PC3", "PC4", 8.71, 5.52) # C_calvus_Atlas - PC3-PC4
p5 <- plot_pca(S_morckhoviensis_Atlas, "PC1", "PC2", 19.07, 13.62) # S_morckhoviensis_Atlas - PC1-PC2
p6 <- plot_pca(S_morckhoviensis_Atlas, "PC3", "PC4", 9.82, 5.93) # S_morckhoviensis_Atlas - PC3-PC4

################################################################################

# Combine plots 

################################################################################

# Combined Plot 
p7 <- plot_grid(
  p1, p2, p3, p4, p5, p6,
  nrow = 2, byrow = FALSE,
  labels = c("A", "B", "C", "D", "E", "F"),
  align="hv")

p7               

##############################################################################

# Measure Euclidian Distances between Original + Alternative Meshes 

##############################################################################

# Create the function 

##############################################################################

# Function to process data, calculate Euclidean distances, and plot
process_data_and_plot <- function(Kernel_data, comparison_data, species_data, species_name, comparison_name, plot_labels) {
  
  # Rename columns in Kernel_data from "PC" to "Comp"
  Kernel_data <- Kernel_data %>%
    rename_with(~ sub("^PC", "Comp", .), starts_with("PC"))
  
  # Merge data into a new dataframe
  merged_data <- merge(Kernel_data, comparison_data, by = "Tip_Label")
  
  # Extract data for 'Arctictis_binturong'
  arctictis_pc_scores <- merged_data %>%
    filter(Tip_Label == "Arctictis_binturong") %>%
    select(starts_with("PC"))
  
  arctictis_comp_scores <- merged_data %>%
    filter(Tip_Label == "Arctictis_binturong") %>%
    select(starts_with("Comp"))
  
  # Calculate Euclidean distances
  data_with_distances <- merged_data %>%
    rowwise() %>%
    mutate(
      PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
      Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
    ) %>%
    ungroup()
  
  # Filter out 'Arctictis_binturong' and select relevant columns
  filtered_data <- data_with_distances %>%
    filter(Tip_Label != "Arctictis_binturong") %>%
    select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)
  
  # Merge with additional species data
  merged_data_filtered <- merge(species_data, filtered_data, by = "Tip_Label")
  
  # Plot data
  p <- ggplot(merged_data_filtered, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
    geom_point(size = 2.5, stroke = 0.4) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
    scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
    scale_color_manual(name = "", values = c("black", "grey60")) +
    scale_fill_manual(values = neworder) +
    labs(
      x = paste0("Arctictis binturong atlas"),
      y = paste0(species_name, " atlas")
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",  # Remove legend
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 11),
      legend.text = element_text(size = 13)
    )
  
  # Perform linear regression
  lm_model <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered)
  lm_summary <- summary(lm_model)
  
  return(list(plot = p, lm_summary = lm_summary))
}

##############################################################################

# Process data and create plots for both comparisons

##############################################################################

# Comparison 1: Arctictis binturong vs Cacajao calvus
result1 <- process_data_and_plot(A_Binturong_Atlas, C_calvus_Atlas, species.data, "Cacajao calvus", "Cacajao calvus", c("A"))
p8 <- result1$plot
lm_model1_summary <- result1$lm_summary

# Comparison 2: Arctictis binturong vs Schizodelphis morckhoviensis
result2 <- process_data_and_plot(A_Binturong_Atlas,S_morckhoviensis_Atlas, species.data, "Schizodelphis morckhoviensis", "Schizodelphis morckhoviensis", c("B"))
p9 <- result2$plot
lm_model2_summary <- result2$lm_summary

# Combined Plot
p10 <- plot_grid(p8, p9, nrow = 1, byrow = FALSE, labels = c("A", "B"), align="hv")
print(p10)

# Print linear regression summaries
print(lm_model1_summary)
print(lm_model2_summary)








