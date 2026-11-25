install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("UpSetR")
install.packages("remotes")
install.packages("ComplexUpset")
remotes::install_url("https://github.com/krassowski/complex-upset/archive/refs/heads/master.zip")

library(ggplot2)
packageVersion("ggplot2")   # should show 3.5.1
library(patchwork)
packageVersion("patchwork")   # 1.2.0
library(ComplexUpset)
packageVersion("ComplexUpset")  # should show 1.3.3
library(RColorBrewer)
library(UpSetR)
library(readxl)

## uPSET PLOT FOR FAECIUM
data_faecium <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/Ecocc_temp.xlsx", sheet = 3)

##Check fata
print(data_faecium)

str(data_faecium)
table(is.na(data_faecium$ERY))
data_faecium[is.na(data_faecium$VAN)==T,]
print(data_faecium)
data <- data_faecium

##ECOFF for faecium TYL is 11, for fecalis is 17
data$ERY = ifelse(data$ERY<= 13.01,1,0)
data$TYL = ifelse(data$TYL<= 11.01,1,0)
data$VAN = ifelse(data$VAN<= 14.01,1,0)
data$CIP = ifelse(data$CIP<= 15.01,1,0)
data$AMP = ifelse(data$AMP<= 16.01,1,0)
data$TET = ifelse(data$TET<= 14.01,1,0)


print(data)

data$ERY[is.na(data$ERY)] <- 0
data$TYL[is.na(data$TYL)] <- 0
data$VAN[is.na(data$VAN)] <- 0
data$CIP[is.na(data$CIP)] <- 0
data$AMP[is.na(data$AMP)] <- 0
data$TET[is.na(data$TET)] <- 0


print(data)

#specifying names used for the antibiotic labels
antibiotics <- c("ERY","TYL", "VAN", "CIP", "AMP", "TET") 
names(antibiotics) <- antibiotics 

print(antibiotics)

#specifying colors to use with the RColorBrewer packagem making them the right length for number of antibiotics 
#could skip this if pallet is continious or has needed number of colours
#colourCount = length(unique(data$FarmID))
#getPalette = colorRampPalette(brewer.pal(3, "BuPu")) #coded down below as specified colours to use
palette(c("#F1C232", "#1F77B4", "#CC0000"))
myPalette <- palette(c("#F1C232", "#1F77B4", "#CC0000"))
myPalette

str(data$FarmID)  # Check the structure of the Farm variable
data$FarmID <- as.factor(data$FarmID)


# constructing the plot with ggplot2 

ComplexUpset::upset(
  data,
  antibiotics,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill=FarmID),
      size = 8,      # Increase dot size in intersection matrix
      linewidth = 4  # Increase line width connecting dots
    ) 
    + geom_bar(stat="identity", width = 0.5)  # colour="black" is replaced with 'width = 1'
    + ylab('Number of Enterococcus faecium isolated')
    + labs(fill = "FarmID") 
    + theme(
      axis.text.y = element_text(size = 10),         # Y-axis text size
      axis.title.x = element_text(size = 10, face = "bold"), # X-axis title size
      axis.title.y = element_text(size = 10, face = "bold"), # Y-axis title size
      legend.text = element_text(size = 10),         # Legend text size
      legend.title = element_text(size = 10)         # Legend title size
    )
  ),
  width_ratio=0.3,
  stripes='white'
) & scale_fill_manual(values = myPalette) 



ggsave("NG_farmUpset_faecium_colors.png", width=11.5, height = 6)
ggsave("NG_farmUpset.svg", width=10, height = 6)

#two way table count
#table(data$AUG, data$TS)

# Create a combination label for each sample
data$Combination <- apply(data[antibiotics], 1, function(x) {
  paste(names(which(x == 1)), collapse = "+")
})
data$Combination[data$Combination == ""] <- "None"

# Summarise which samples belong to each combination
samples_per_combination <- aggregate(IsolateID ~ Combination, data = data, FUN = paste, collapse = ", ")

# Count how many samples in each combination
samples_per_combination$Count <- sapply(strsplit(samples_per_combination$IsolateID, ", "), length)

# View and save the table
print(samples_per_combination)
write.csv(samples_per_combination, "samples_per_combination_faecium.csv", row.names = FALSE)

table(data$FarmID, useNA = "always")


samples_per_combination <- aggregate(IsolateID ~ Combination, data = data, FUN = paste, collapse = ", ")
samples_per_combination$Count <- sapply(strsplit(samples_per_combination$IsolateID, ", "), length)
print(samples_per_combination)

install.packages("writexl")
library(writexl)

write_xlsx(samples_per_combination, "samples_per_combination_faecium.xlsx")
