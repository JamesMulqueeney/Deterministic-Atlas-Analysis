# Paper - Assessing the application of landmark-free morphometrics to macroevolutionary analyses

# Author: James M. Mulqueeney 

# Date Last Modified: 28/02/2025

# Plot Evolutionary Rates & Disparity Data Results 

#########################################################################################

# Load required libraries
library(ggplot2)

# Define input directory
input_dir <- "C:/Users/jmm1e21/OneDrive - University of Southampton/Documents/Main PhD Work/Chapter 3 - Comparison of Automated Methods/Final/BMC Ecology and Evolution/Data & Code/Data/Data/"

#########################################################################################

# Disparity Plots 

#########################################################################################

# Diet  

# Read in the data 
diet <- read.csv(file.path(input_dir, "Data_A18-Diet_Disparity_Rates_v1.csv"))

D1 <- ggplot(diet, aes(x = log(M_Disparity), y = log(K_Disparity) , fill = Type, colour = Type)) +
  geom_point(size = 4, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, size = 1,  aes(fill = Type)) +
  labs(x = expression(paste("log"[10], " (Manual Landmarking Disparity)")),
    y = expression(paste("log"[10], " (DAA Disparity)")), color = "Kernel Width") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"))
D1  

# Group by Type and calculate Spearman's rank correlation for each group
spearman_results <- diet %>%
  group_by(Type) %>%
  summarise(
    spearman_corr = cor(log(M_Disparity), log(K_Disparity), method = "spearman", use = "complete.obs"),
    p_value = cor.test(log(M_Disparity), log(K_Disparity), method = "spearman")$p.value  # Extract p-value
  )

# View the Spearman's correlation and p-values
print(spearman_results)

#########################################################################################

# Locomotion

# Read in the data 
loc <- read.csv(file.path(input_dir, "Data_A19-Locomotion_Disparity_Rates_v1.csv"))

D2 <- ggplot(loc, aes(x = log(M_Disparity), y = log(K_Disparity), fill = Type, colour = Type)) +
  geom_point(size = 4, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, size = 1, aes(fill = Type)) +
  labs(x = expression(paste("log"[10], " (Manual Landmarking Disparity)")),
    y = expression(paste("log"[10], " (DAA Disparity)")), color = "Kernel Width") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"))
D2 

# Group by Type and calculate Spearman's rank correlation for each group
spearman_results2 <- loc %>%
  group_by(Type) %>%
  summarise(
    spearman_corr = cor(log(M_Disparity), log(K_Disparity), method = "spearman", use = "complete.obs"),
    p_value = cor.test(log(M_Disparity), log(K_Disparity), method = "spearman")$p.value  # Extract p-value
  )

# View the Spearman's correlation and p-values
print(spearman_results2)

#########################################################################################

# Rates Plots 

#########################################################################################

# Diet  

R1 <- ggplot(diet, aes(x = log(M_Rate), y = log(K_Rate) , fill = Type, colour = Type)) +
  geom_point(size = 4, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, size = 1, aes(fill = Type)) +
  labs(x = expression(paste("log"[10], " (Manual Landmarking Rate)")),
    y = expression(paste("log"[10], " (DAA Rate)")), color = "Kernel Width") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"))
R1 

# Group by Type and calculate Spearman's rank correlation for each group
spearman_results1 <- diet %>%
  group_by(Type) %>%
  summarise(
    spearman_corr = cor(log(M_Rate), log(K_Rate), method = "spearman", use = "complete.obs"),
    p_value = cor.test(log(M_Rate), log(K_Rate), method = "spearman")$p.value  # Extract p-value
  )

print(spearman_results1)

#########################################################################################

# Locomotion  

R2 <- ggplot(loc, aes(x = log(M_Rate), y = log(K_Rate) , fill = Type, colour = Type)) +
  geom_point(size = 4, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, size = 1, aes(fill = Type)) +
  labs(x = expression(paste("log"[10], " (Manual Landmarking Rate)")),
    y = expression(paste("log"[10], " (DAA Rate)")), color = "Kernel Width") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"))
R2 

# Group by Type and calculate Spearman's rank correlation for each group
spearman_results2 <- loc %>%
  group_by(Type) %>%
  summarise(
    spearman_corr = cor(log(M_Rate), log(K_Rate), method = "spearman", use = "complete.obs"),
    p_value = cor.test(log(M_Rate), log(K_Rate), method = "spearman")$p.value  # Extract p-value
  )

print(spearman_results2)

#########################################################################################

# Combined Plot 

#########################################################################################

g1 <- plot_grid(D1, D2, R1, R2,  
                nrow = 2, byrow = FALSE,
                labels = c("A", "B", "C", "D"),
                align="hv")

g1
