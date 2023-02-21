# Installs -----------------------------------------------------------------
install.packages("devtools")
install.packages(c("dplyr", "tidyr", "ggplot2", "MASS", "meta", "ggrepel", "DT"))
devtools::install_github("PheWAS/PheWAS")
library(PheWAS)

# Dataloading -------------------------------------------------------------
setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")
data <- read.csv("data/phewas_tables/TRPV1_phewas_table.csv")
phenotype <- data[, c(5:ncol(data))]
phenotype <- createPhenotypes(phenotype)

