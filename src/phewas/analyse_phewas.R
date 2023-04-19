library(data.table)
library(tidyverse)
library(dplyr)
library(PheWAS)

setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")

phewas_list <- readRDS("outputs/phe_results_v1_2.rds")

combined_phewas <- bind_rows(phewas_list)
combined_phewas <- setDT(combined_phewas)
combined_phewas <- combined_phewas[combined_phewas$p < 0.05,]
combined_phewas <- addPhecodeInfo(combined_phewas)
