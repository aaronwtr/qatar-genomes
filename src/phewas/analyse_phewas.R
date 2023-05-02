library(data.table)
library(tidyverse)
library(dplyr)
library(PheWAS)

setwd("~/Desktop/PhD/Research/QMUL/Research/qatar-genomes")

phewas_list <- readRDS("outputs/full_results/phewas_age_gender_adjusted.rds")
phewas_list_bmi <- readRDS("outputs/full_results/phewas_age_gender_bmi_adjusted.rds")

combined_phewas <- bind_rows(phewas_list)
combined_phewas <- setDT(combined_phewas)
combined_phewas <- combined_phewas[combined_phewas$p < 0.05,]
combined_phewas <- addPhecodeInfo(combined_phewas)

combined_phewas_bmi <- bind_rows(phewas_list_bmi)
combined_phewas_bmi <- setDT(combined_phewas_bmi)
combined_phewas_bmi <- combined_phewas_bmi[combined_phewas_bmi$p < 0.05,]
combined_phewas_bmi <- addPhecodeInfo(combined_phewas_bmi)
