library(data.table)
library(tidyverse)
library(dplyr)
library(PheWAS)

## Load and preprocess data
# TODO Calculate bonferoni adjusted p values in put in separate column
setwd("~/Desktop/PhD/Research/QMUL/Research/qatar-genomes")

phewas_list <- readRDS("outputs/full_results/phewas_age_gender_adjusted.rds")
phewas_list_bmi <- readRDS("outputs/full_results/phewas_age_gender_bmi_adjusted.rds")

combined_phewas <- bind_rows(phewas_list)
combined_phewas <- setDT(combined_phewas)
combined_phewas <- combined_phewas[combined_phewas$p < 0.05,]
combined_phewas <- addPhecodeInfo(combined_phewas)
group_factor <- factor(combined_phewas$group)
group_num <- as.numeric(group_factor)
combined_phewas$groupnum <- group_num


combined_phewas_bmi <- bind_rows(phewas_list_bmi)
combined_phewas_bmi <- setDT(combined_phewas_bmi)
combined_phewas_bmi <- combined_phewas_bmi[combined_phewas_bmi$p < 0.05,]
combined_phewas_bmi <- addPhecodeInfo(combined_phewas_bmi)
group_factor <- factor(combined_phewas_bmi$group)
group_num <- as.numeric(group_factor)
combined_phewas_bmi$groupnum <- group_num

## Find out the most occurring, significant variants 
snp_freq <- table(combined_phewas$snp)
snp_freq_sorted <- data.table(sort(snp_freq, decreasing = TRUE))
setnames(snp_freq_sorted, "V1", "snp")
ith_snp <- snp_freq_sorted[10, "snp"]


significant_hits_variant_i <- combined_phewas[snp == as.character(ith_snp),]
pheno_desc <- significant_hits_variant_i[, .(phenotype, description)]

pdf_file <- sprintf("plots/phewas_plots/manhattan_plots/%s.pdf", ith_snp)

# Open a new PDF device
pdf(file = pdf_file)

# Generate the phewas plot
phewasManhattan(d=significant_hits_variant_i, 
                annotate.phenotype.description=pheno_desc, x.group.labels=F,
                title=ith_snp, annotate.angle=0, annotate.size=5)

dev.off()