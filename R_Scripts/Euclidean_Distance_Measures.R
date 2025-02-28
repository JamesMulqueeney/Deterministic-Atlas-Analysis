# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Euclidean Distance Comparisons 

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

# Extract data for 'Arctictis_binturong'
arctictis_pc_scores <- merged_data1 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("PC"))

arctictis_comp_scores <- merged_data1 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("Comp"))

# Calculate Euclidean distances for all PC and Comp components
data_with_distances1 <- merged_data1 %>%
  rowwise() %>%
  mutate(
    PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
    Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
  ) %>%
  ungroup()

# Filter out 'Arctictis_binturong' and select relevant columns
filtered_data1 <- data_with_distances1 %>%
  filter(Tip_Label != "Arctictis_binturong") %>%
  select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)

# Merge with additional species data (assuming 'species_data' is correctly defined)
merged_data_filtered1 <- merge(species.data, filtered_data1, by = "Tip_Label")

# Plotting the data
p1 <- ggplot(merged_data_filtered1, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = "Manual Landmarking",
    y = "Deterministic Atlas Analysis"
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

print(p1)

# Perform linear regression
lm_model1 <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered1)

# Extract summary statistics
summary(lm_model1)

# Filter out Cetacea from the dataset
merged_data_no_cetacea1 <- merged_data_filtered1 %>%
  filter(Order != "Cetacea")

# Perform linear regression excluding cetaceans
lm_model_no_cetacea1 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea1)

# Extract summary statistics
summary(lm_model_no_cetacea1)

# Filter out Primates from the dataset
merged_data_no_primates1 <- merged_data_filtered1 %>%
  filter(Order != "Primates")

# Perform linear regression excluding cetaceans
lm_model_no_primates1 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_primates1)

# Extract summary statistics
summary(lm_model_no_primates1)

# Filter out Cetacea and Primates from the dataset
merged_data_no_cetacea_primates1 <- merged_data_filtered1 %>%
  filter(Order != "Cetacea", Order != "Primates")

# Perform linear regression excluding Cetacea and Primates
lm_model_no_cetacea_primates1 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea_primates1)

# Extract summary statistics
summary(lm_model_no_cetacea_primates1)

## Extract for Orders 

# Filter for Orders with at least 10 specimens
merged_data_filtered1 <- merged_data_filtered1 %>%
  group_by(Order) %>%
  filter(n() >= 10) %>%
  ungroup()

# Group by Superorder and fit the linear model for each group
lm_results1 <- merged_data_filtered1 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance, data = .))

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value1 <- lm_results1 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"]  # Extract p-value for Comp_distance
  )

# View the R-squared and p-values
print(lm_r_squared_p_value1)

################################################################################

# Aligned_Only: Manual vs Kernel 20.0 (Control points = 270 )

################################################################################

# Merge the data into a new dataframe 
merged_data2 <- merge(pcscores, Kernel_20.0, by = "Tip_Label")

# Extract data for 'Arctictis_binturong'
arctictis_pc_scores <- merged_data2 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("PC"))

arctictis_comp_scores <- merged_data2 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("Comp"))

# Calculate Euclidean distances for all PC and Comp components
data_with_distances2 <- merged_data2 %>%
  rowwise() %>%
  mutate(
    PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
    Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
  ) %>%
  ungroup()

# Filter out 'Arctictis_binturong' and select relevant columns
filtered_data2 <- data_with_distances2 %>%
  filter(Tip_Label != "Arctictis_binturong") %>%
  select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)

# Merge with additional species data (assuming 'species_data' is correctly defined)
merged_data_filtered2 <- merge(species.data, filtered_data2, by = "Tip_Label")

# Plotting the data
p2 <- ggplot(merged_data_filtered2, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = "Manual Landmarking",
    y = "Deterministic Atlas Analysis"
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

print(p2)

# Perform linear regression
lm_model2 <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered2)

# Extract summary statistics
summary(lm_model2)

# Filter out Cetacea from the dataset
merged_data_no_cetacea2 <- merged_data_filtered2 %>%
  filter(Order != "Cetacea")

# Perform linear regression excluding cetaceans
lm_model_no_cetacea2 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea2)

# Extract summary statistics
summary(lm_model_no_cetacea2)

# Filter out Primates from the dataset
merged_data_no_primates2 <- merged_data_filtered2 %>%
  filter(Order != "Primates")

# Perform linear regression excluding cetaceans
lm_model_no_primates2 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_primates2)

# Extract summary statistics
summary(lm_model_no_primates2)

# Filter out Cetacea and Primates from the dataset
merged_data_no_cetacea_primates2 <- merged_data_filtered2 %>%
  filter(Order != "Cetacea", Order != "Primates")

# Perform linear regression excluding Cetacea and Primates
lm_model_no_cetacea_primates2 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea_primates2)

# Extract summary statistics
summary(lm_model_no_cetacea_primates2)

# Extract for Orders 

# Filter for Orders with at least 10 specimens
merged_data_filtered2 <- merged_data_filtered2 %>%
  group_by(Order) %>%
  filter(n() >= 10) %>%
  ungroup()

# Group by Superorder and fit the linear model for each group
lm_results2 <- merged_data_filtered2 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance, data = .))

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value2 <- lm_results2 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"]  # Extract p-value for Comp_distance
  )

# View the R-squared and p-values
print(lm_r_squared_p_value2)

################################################################################

# Aligned_Only: Manual vs Kernel 10.0 (Control points = 1782 )

################################################################################

# Merge the data into a new dataframe 
merged_data3 <- merge(pcscores, Kernel_10.0, by = "Tip_Label")

# Extract data for 'Arctictis_binturong'
arctictis_pc_scores <- merged_data3 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("PC"))

arctictis_comp_scores <- merged_data3 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("Comp"))

# Calculate Euclidean distances for all PC and Comp components
data_with_distances3 <- merged_data3 %>%
  rowwise() %>%
  mutate(
    PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
    Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
  ) %>%
  ungroup()

# Filter out 'Arctictis_binturong' and select relevant columns
filtered_data3 <- data_with_distances3 %>%
  filter(Tip_Label != "Arctictis_binturong") %>%
  select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)

# Merge with additional species data (assuming 'species_data' is correctly defined)
merged_data_filtered3 <- merge(species.data, filtered_data3, by = "Tip_Label")

# Plotting the data
p3 <- ggplot(merged_data_filtered3, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = "Manual Landmarking",
    y = "Deterministic Atlas Analysis"
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

print(p3)

# Perform linear regression
lm_model3 <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered3)

# Extract summary statistics
summary(lm_model3)

# Filter out Cetacea from the dataset
merged_data_no_cetacea3 <- merged_data_filtered3 %>%
  filter(Order != "Cetacea")

# Perform linear regression excluding cetaceans
lm_model_no_cetacea3 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea3)

# Extract summary statistics
summary(lm_model_no_cetacea3)

# Filter out Primates from the dataset
merged_data_no_primates3 <- merged_data_filtered3 %>%
  filter(Order != "Primates")

# Perform linear regression excluding cetaceans
lm_model_no_primates3 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_primates3)

# Extract summary statistics
summary(lm_model_no_primates3)

# Filter out Cetacea and Primates from the dataset
merged_data_no_cetacea_primates3 <- merged_data_filtered3 %>%
  filter(Order != "Cetacea", Order != "Primates")

# Perform linear regression excluding Cetacea and Primates
lm_model_no_cetacea_primates3 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea_primates3)

# Extract summary statistics
summary(lm_model_no_cetacea_primates3)

# Extract for Orders 

# Filter for Orders with at least 10 specimens
merged_data_filtered3 <- merged_data_filtered3 %>%
  group_by(Order) %>%
  filter(n() >= 10) %>%
  ungroup()

# Group by Superorder and fit the linear model for each group
lm_results3 <- merged_data_filtered3 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance, data = .))

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value3 <- lm_results3 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"]  # Extract p-value for Comp_distance
  )

# View the R-squared and p-values
print(lm_r_squared_p_value3)

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

# Poisson Mesh: Manual vs Kernel 40.0 (Control points = 45)

################################################################################

# Merge the data into a new dataframe 
merged_data4 <- merge(pcscores, Kernel_40.0_P, by = "Tip_Label")

# Extract data for 'Arctictis_binturong'
arctictis_pc_scores <- merged_data4 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("PC"))

arctictis_comp_scores <- merged_data4 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("Comp"))

# Calculate Euclidean distances for all PC and Comp components
data_with_distances4 <- merged_data4 %>%
  rowwise() %>%
  mutate(
    PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
    Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
  ) %>%
  ungroup()

# Filter out 'Arctictis_binturong' and select relevant columns
filtered_data4 <- data_with_distances4 %>%
  filter(Tip_Label != "Arctictis_binturong") %>%
  select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)

# Merge with additional species data (assuming 'species_data' is correctly defined)
merged_data_filtered4 <- merge(species.data, filtered_data4, by = "Tip_Label")

# Plotting the data
p4 <- ggplot(merged_data_filtered4, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = "Manual Landmarking",
    y = "Deterministic Atlas Analysis"
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

print(p4)

# Perform linear regression
lm_model4 <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered4)

# Extract summary statistics
summary(lm_model4)

# Filter out Cetacea from the dataset
merged_data_no_cetacea4 <- merged_data_filtered4 %>%
  filter(Order != "Cetacea")

# Perform linear regression excluding cetaceans
lm_model_no_cetacea4 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea4)

# Extract summary statistics
summary(lm_model_no_cetacea4)

# Filter out Primates from the dataset
merged_data_no_primates4 <- merged_data_filtered4 %>%
  filter(Order != "Primates")

# Perform linear regression excluding cetaceans
lm_model_no_primates4 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_primates4)

# Extract summary statistics
summary(lm_model_no_primates4)

# Filter out Cetacea and Primates from the dataset
merged_data_no_cetacea_primates4 <- merged_data_filtered4 %>%
  filter(Order != "Cetacea", Order != "Primates")

# Perform linear regression excluding Cetacea and Primates
lm_model_no_cetacea_primates4 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea_primates4)

# Extract summary statistics
summary(lm_model_no_cetacea_primates4)

# Extract for Orders 

# Filter for Orders with at least 10 specimens
merged_data_filtered4  <- merged_data_filtered4 %>%
  group_by(Order) %>%
  filter(n() >= 10) %>%
  ungroup()

# Group by Superorder and fit the linear model for each group
lm_results4 <- merged_data_filtered4 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance, data = .))

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value4 <- lm_results4 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"]  # Extract p-value for Comp_distance
  )

# View the R-squared and p-values
print(lm_r_squared_p_value4)

# Group by Order and fit the linear model for each group including Status as a predictor
lm_results4 <- merged_data_filtered4 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance + Status, data = .))  # Add Status as a predictor

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value4 <- lm_results4 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"],  # Extract p-value for Comp_distance
    p_value_status = summary(lm_fit)$coefficients["StatusExtinct", "Pr(>|t|)"]  # Extract p-value for Status
  )

# View the R-squared and p-values
print(lm_r_squared_p_value4)

################################################################################

# Poisson Mesh: Manual vs Kernel 20.0 (Control points = 270)

################################################################################

# Merge the data into a new dataframe 
merged_data5 <- merge(pcscores, Kernel_20.0_P, by = "Tip_Label")

# Extract data for 'Arctictis_binturong'
arctictis_pc_scores <- merged_data5 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("PC"))

arctictis_comp_scores <- merged_data5 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("Comp"))

# Calculate Euclidean distances for all PC and Comp components
data_with_distances5 <- merged_data5 %>%
  rowwise() %>%
  mutate(
    PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
    Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
  ) %>%
  ungroup()

# Filter out 'Arctictis_binturong' and select relevant columns
filtered_data5 <- data_with_distances5 %>%
  filter(Tip_Label != "Arctictis_binturong") %>%
  select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)

# Merge with additional species data (assuming 'species_data' is correctly defined)
merged_data_filtered5 <- merge(species.data, filtered_data5, by = "Tip_Label")

# Plotting the data
p5 <- ggplot(merged_data_filtered5, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = "Manual Landmarking",
    y = "Deterministic Atlas Analysis"
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

print(p5)

# Perform linear regression
lm_model5 <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered5)

# Extract summary statistics
summary(lm_model5)

# Filter out Cetacea from the dataset
merged_data_no_cetacea5 <- merged_data_filtered5 %>%
  filter(Order != "Cetacea")

# Perform linear regression excluding cetaceans
lm_model_no_cetacea5 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea5)

# Extract summary statistics
summary(lm_model_no_cetacea5)

# Filter out Primates from the dataset
merged_data_no_primates5 <- merged_data_filtered5 %>%
  filter(Order != "Primates")

# Perform linear regression excluding cetaceans
lm_model_no_primates5 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_primates5)

# Extract summary statistics
summary(lm_model_no_primates5)

# Filter out Cetacea and Primates from the dataset
merged_data_no_cetacea_primates5 <- merged_data_filtered5 %>%
  filter(Order != "Cetacea", Order != "Primates")

# Perform linear regression excluding Cetacea and Primates
lm_model_no_cetacea_primates5 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea_primates5)

# Extract summary statistics
summary(lm_model_no_cetacea_primates5)

# Extract for Orders 

# Filter for Orders with at least 10 specimens
merged_data_filtered5 <- merged_data_filtered5 %>%
  group_by(Order) %>%
  filter(n() >= 10) %>%
  ungroup()

# Group by Superorder and fit the linear model for each group
lm_results5 <- merged_data_filtered5 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance, data = .))

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value5 <- lm_results5 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"]  # Extract p-value for Comp_distance
  )

# View the R-squared and p-values
print(lm_r_squared_p_value5)

################################################################################

# Poisson Mesh: Manual vs Kernel 10.0 (Control points = 1782)

################################################################################

# Merge the data into a new dataframe 
merged_data6 <- merge(pcscores, Kernel_10.0_P, by = "Tip_Label")

# Extract data for 'Arctictis_binturong'
arctictis_pc_scores <- merged_data6 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("PC"))

arctictis_comp_scores <- merged_data6 %>%
  filter(Tip_Label == "Arctictis_binturong") %>%
  select(starts_with("Comp"))

# Calculate Euclidean distances for all PC and Comp components
data_with_distances6 <- merged_data6 %>%
  rowwise() %>%
  mutate(
    PC_distance = sqrt(sum((c_across(starts_with("PC")) - arctictis_pc_scores[1, ])^2)),
    Comp_distance = sqrt(sum((c_across(starts_with("Comp")) - arctictis_comp_scores[1, ])^2))
  ) %>%
  ungroup()

# Filter out 'Arctictis_binturong' and select relevant columns
filtered_data6 <- data_with_distances6 %>%
  filter(Tip_Label != "Arctictis_binturong") %>%
  select(Tip_Label, starts_with("PC"), starts_with("Comp"), PC_distance, Comp_distance)

# Merge with additional species data (assuming 'species_data' is correctly defined)
merged_data_filtered6 <- merge(species.data, filtered_data6, by = "Tip_Label")

# Plotting the data
p6 <- ggplot(merged_data_filtered6, aes(x = Comp_distance, y = PC_distance, shape = Mesh.Type, fill = Order, color = Status)) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "red", size = 1) +  
  scale_shape_manual(name = "", values = c(21, 22, 23, 25, 24)) +
  scale_color_manual(name = "", values = c("black", "grey60")) +
  scale_fill_manual(values = neworder) +
  labs(
    x = "Manual Landmarking",
    y = "Deterministic Atlas Analysis"
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

print(p6)

# Perform linear regression
lm_model6 <- lm(PC_distance ~ Comp_distance, data = merged_data_filtered6)

# Extract summary statistics
summary(lm_model6)

# Filter out Cetacea from the dataset
merged_data_no_cetacea6 <- merged_data_filtered6 %>%
  filter(Order != "Cetacea")

# Perform linear regression excluding cetaceans
lm_model_no_cetacea6 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea6)

# Extract summary statistics
summary(lm_model_no_cetacea6)

# Filter out Cetacea from the dataset
merged_data_no_primates6 <- merged_data_filtered6 %>%
  filter(Order != "Primates")

# Perform linear regression excluding cetaceans
lm_model_no_primates6 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_primates6)

# Extract summary statistics
summary(lm_model_no_primates6)

# Filter out Cetacea and Primates from the dataset
merged_data_no_cetacea_primates6 <- merged_data_filtered6 %>%
  filter(Order != "Cetacea", Order != "Primates")

# Perform linear regression excluding Cetacea and Primates
lm_model_no_cetacea_primates6 <- lm(PC_distance ~ Comp_distance, data = merged_data_no_cetacea_primates6)

# Extract summary statistics
summary(lm_model_no_cetacea_primates6)

# Filter for Orders with at least 10 specimens
merged_data_filtered6  <- merged_data_filtered6 %>%
  group_by(Order) %>%
  filter(n() >= 10) %>%
  ungroup()

# Group by Superorder and fit the linear model for each group
lm_results6 <- merged_data_filtered6 %>%
  group_by(Order) %>%
  do(lm_fit = lm(PC_distance ~ Comp_distance, data = .))

# Calculate and extract R-squared and p-values for each group
lm_r_squared_p_value6 <- lm_results6 %>%
  summarise(
    Order, 
    r_squared = summary(lm_fit)$r.squared,
    p_value = summary(lm_fit)$coefficients["Comp_distance", "Pr(>|t|)"]  # Extract p-value for Comp_distance
  )

# View the R-squared and p-values
print(lm_r_squared_p_value6)

################################################################################

# Combined plot 

################################################################################

# Combined Plot 
c1 <- plot_grid(
  p1,p2,p3,p4,p5,p6,
  nrow = 2, byrow = TRUE,
  labels = c("A", "B", "C", "D", "E", "F"),
  align="hv")

c1

################################################################################

# T-Test Comparison of Aligned-only vs Poisson Results

################################################################################

# Define Euclidean statistic values
aligned <- c(0.1469, 0.1262, 0.1714) # Mean =  0.285 SD = 0.0229
poisson <- c(0.5317, 0.4343, 0.4245) # Mean =  0.619 SD = 0.0354

# Measure Means  
signif(mean(aligned), 3)
signif(sd(aligned), 3)

# Measure Standard Deviation 
signif(mean(poisson), 3)
signif(sd(poisson), 3)

# Perform paired t-test
t.test(poisson, aligned, paired = TRUE)