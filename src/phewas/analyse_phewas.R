library(data.table)
library(tidyverse)
library(dplyr)
library(PheWAS)
library(ggplot2)
library(writexl)

## Load and preprocess data
# TODO Calculate bonferoni adjusted p values in put in separate column
setwd("~/Desktop/PhD/Research/QMUL/Research/qatar-genomes")

phewas_list <- readRDS("outputs/full_results/phewas_age_gender_adjusted.rds")
phewas_list_bmi <- readRDS("outputs/full_results/phewas_age_gender_bmi_adjusted.rds")

combined_phewas <- bind_rows(phewas_list)
combined_phewas <- setDT(combined_phewas)
combined_phewas <- combined_phewas[combined_phewas$p < 0.05,]
combined_phewas <- addPhecodeInfo(combined_phewas, groupnums=T)

## Generating list of results, lowest association for each variant

phewas_sing_snp <- combined_phewas %>%
  group_by(snp) %>%
  summarize(p = min(p), description = first(description), 
            phenotype = first(phenotype)) %>%
  ungroup()

write_xlsx(phewas_sing_snp, 
           "outputs/full_results/phewas_summary_age_gender_adjusted.xlsx")

## Plotting ouput phecode group distribution

group_counts <- count(combined_phewas, group)
group_counts <- group_counts[order(-group_counts$n), ]

chart <- ggplot(group_counts, aes(x = reorder(group, -n), y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL) 

print(chart)

## Load and preprecess BMI output

combined_phewas_bmi <- bind_rows(phewas_list_bmi)
combined_phewas_bmi <- setDT(combined_phewas_bmi)
combined_phewas_bmi <- combined_phewas_bmi[combined_phewas_bmi$p < 0.05,]
combined_phewas_bmi <- addPhecodeInfo(combined_phewas_bmi, groupnums=T)

## Find out the most occurring, significant variants 
snp_freq <- table(combined_phewas$snp)
snp_freq_sorted <- data.table(sort(snp_freq, decreasing = TRUE))
setnames(snp_freq_sorted, "V1", "snp")
# ith_snp <- snp_freq_sorted[1, "snp"]
ith_snp <- "GLB1L3_11.134162108.G.A"

significant_hits_variant_i <- combined_phewas[snp == as.character(ith_snp),]
pheno_desc <- significant_hits_variant_i[, .(phenotype, description)]

pdf_file <- sprintf("plots/phewas_plots/manhattan_plots/%s.pdf", ith_snp)

# Open a new PDF device
pdf(file = pdf_file)

# Generate the phewas plot
phewasManhattan(d=significant_hits_variant_i, 
                annotate.phenotype.description=pheno_desc, x.group.labels=F,
                title=ith_snp, annotate.angle=0, annotate.size=3)

dev.off()