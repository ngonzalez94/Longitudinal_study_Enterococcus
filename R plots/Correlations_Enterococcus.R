library(hexbin)
#library(RColorBrewer)
library(tidyverse)
# library(viridis)
library(GGally)
library(ggplot2)

library(ggcorrplot)
#library(ggstatsplot)


#packageVersion("tidyverse")
#packageVersion("GGally")
#packageVersion("ggplot2")
#packageVersion("RColorBrewer")
#packageVersion("hexbin")

## setwd()
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/")

library(readxl)

# Assuming Correlation.xlsx is in the current working directory and you want to read the first sheet (#1 faecalis # 3faecium)
Correlations <- read_excel("Ecocc_temp.xlsx", sheet = 1)

## load in data
#Correlation <- read.table("Correlations.xlsx", header = TRUE)

print(Correlations)
Correlations <- Correlations[-1]

# Check the structure and summary of the data
print(str(Correlations))
print(summary(Correlations))

#Remove NA 
# Remove rows with NA values
Correlations_clean <- na.omit(Correlations)

## run function
ggpairs_hex <- function(df, hexbins = 30) {
  # REF: https://stackoverflow.com/questions/20872133/using-stat-binhex-with-ggpairs
  p <- ggpairs(df, lower="blank",
               upper = list(continuous = wrap("cor", method= "spearman", color = "black")))
  seq <- 1:ncol(df)
  for (x in seq)
    for (y in seq) 
      if (y>x) 
        p <- putPlot(p, 
                     ggplot(df, aes_string(x=names(df)[x], y=names(df)[y])) + 
                       stat_binhex(bins = hexbins) + 
                       scale_fill_gradient(low='coral', high='darkblue'),
                     #scale_fill_brewer(palette="PuOr"),
                     #scale_fill_gradient(colours = my_colors(11)), 
                     y,  x)
  return(p)
}


ggpairs_hex(Correlations_clean[,1:6]) 
ggtitle("Correlation between the diffusion zone's sizes of antimicrobails tested against Enterococcus (E. faecalis")
ggsave("Correlations_All.pdf", width = 28, height = 19, units = "cm")

#############################

# correlation plots(ignore this)

class(BVs_FPrime)
samp <- BVs_FPrime %>% remove_rownames %>% column_to_rownames(var="Animal")
ggcorrplot(cor(samp))
ggsave("BVs_FPrime_corrplot.pdf", width = 15, height = 12, units = "cm")

