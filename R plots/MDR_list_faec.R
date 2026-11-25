# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(rlang)
#install.packages("writexl")
library(writexl)
library(openxlsx)


## setwd()
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/")

library(readxl)

# Assuming Correlation.xlsx is in the current working directory and you want to read the first sheet (sheet 5 faecium,, sheet 6 faecalis)
data <- read_excel("Ecocc_temp.xlsx", sheet = 5)

# Define the cut-off values for each antibiotic (replace with your actual cut-off values)(ECOFF FAECIUM = 11)(ECOFF FAECALIS = 17)
cut_offs <- c(
  AMP = 16.01,
  CIP = 15.01,
  ERY = 13.01,
  TET = 14.01,
  TYL = 11.01,
  VAN = 14.01
  )

# Create a new column to count the number of resistances
data <- data %>%
  rowwise() %>%
  mutate(ResistanceCount = sum(c_across(names(cut_offs)) < unlist(cut_offs)))

# Filter isolates that are resistant to 3 or more antibiotics
resistant_3_or_more <- data %>%
  filter(ResistanceCount >= 3)

# Select the columns for output (Isolate, Farm, and antibiotic resistance zones)
final_output <- resistant_3_or_more %>%
  select(IsolateID, FarmID, everything())

# Print the result
print(final_output)

View(final_output)

#Save the result into a new Excel file
output_file_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/resistant_isolates_highlighted_faecIUM.xlsx"
write_xlsx(final_output,output_file_path)

