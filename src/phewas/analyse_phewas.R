library(data.table)
library(tidyverse)
library(dplyr)
library(PheWAS)

setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")

phewas_list <- readRDS("outputs/phe_results_v1.rds")
trpv1 <- phewas_list[["TRPV1_17.3493179.C.CGG"]]
results_phenotest_info <- addPhecodeInfo(trpv1)

combined_phewas <- bind_rows(phewas_list)
combined_phewas <- setDT(combined_phewas)
combined_phewas <- combined_phewas[combined_phewas$p < 0.05,]

