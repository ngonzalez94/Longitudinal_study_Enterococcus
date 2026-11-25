# Load necessary libraries
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

# Read data from Excel file
data <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/Enterococcus_FF.xlsx", sheet = 6)

# Reshape data for plotting
data_long <- pivot_longer(data, cols = -SampleID, names_to = "Antimicrobial", values_to = "Diffusion_Zone")

#Calculate ECOFF as median diffusion zone for each antimicrobial
ecoff <- data_long %>%
  group_by(Antimicrobial) %>%
  summarise(ECOFF = median(Diffusion_Zone, na.rm = TRUE)) 

# Plot ECOFFs
ggplot(data_long, aes(x = Diffusion_Zone)) +
  geom_density() + 
  facet_wrap(~Antimicrobial, scales = "free") +
  geom_vline(data = ecoff, aes(xintercept = ECOFF, color = Antimicrobial), linetype = "dashed") +
  labs(title = "Distribution of Diffusion Zones by Antimicrobial Faecium",
       x = "Diffusion Zone",
       y = "Density") +
  theme_minimal()

data <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/Enterococcus_FF.xlsx", sheet = 4)

# Reshape data for plotting
data_long <- pivot_longer(data, cols = -SampleID, names_to = "Antimicrobial", values_to = "Diffusion_Zone")

#Calculate ECOFF as median diffusion zone for each antimicrobial
ecoff <- data_long %>%
  group_by(Antimicrobial) %>%
  summarise(ECOFF = median(Diffusion_Zone, na.rm = TRUE)) 

# Plot ECOFFs
ggplot(data_long, aes(x = Diffusion_Zone)) +
  geom_density() + 
  facet_wrap(~Antimicrobial, scales = "free") +
  geom_vline(data = ecoff, aes(xintercept = ECOFF, color = Antimicrobial), linetype = "dashed") +
  labs(title = "Distribution of Diffusion Zones by Antimicrobial Faecalis",
       x = "Diffusion Zone",
       y = "Density") +
  theme_minimal()
