# Installs -----------------------------------------------------------------
install.packages("devtools")
install.packages(c("dplyr", "tidyr", "ggplot2", "MASS", "meta", "ggrepel", "DT"))
install.packages("readxl")
devtools::install_github("PheWAS/PheWAS")
library(PheWAS)
library(readxl)


# Dataloading -------------------------------------------------------------
setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")
data <- read_excel("data/phewas_tables/TRPV1_phewas_table.xlsx", sheet = 2)
phenotype <- data[, c(5:ncol(data))]
phenotype <- createPhenotypes(phenotype)
