library(PheWAS)
library(data.table)
library(tidyverse)
library(readxl)
library(dplyr)
library(progress)

setwd("~/Desktop/PhD/Research/QMUL/Research/qatar-genomes")

genodata <- data.frame(fread("data/phewas_input/gene_data_bmi.csv"))
phenodata <- fread("data/phewas_input/phenotype_data.csv")

## Processing phenotype data

id.sex <- genodata[,c("patient_id","gender")]
id.sex$gender[which(id.sex$gender=="MALE")]<-"M"
id.sex$gender[which(id.sex$gender=="FEMALE")]<-"F"

genodata$patient_id <- as.character(genodata$patient_id)
genodata$gender <- as.factor(genodata$gender)
genodata$BMI <- as.numeric(genodata$BMI)
id.sex$patient_id <- as.character(id.sex$patient_id)
id.sex$gender <- as.factor(id.sex$gender)

phenodata <- setDT(phenodata)
phenodata<- phenodata[, patient_id := as.character(patient_id)]
id.sex.test <- filter(id.sex, patient_id %in% phenodata$patient_id)

phenotypes <- createPhenotypes(phenodata, min.code.count = 1, 
                               id.sex = id.sex.test, 
                               vocabulary.map = PheWAS::phecode_map_icd10,
                               rollup.map=PheWAS::phecode_rollup_map)
phenotypes_nonan<- phenotypes
phenotypes_nonan[is.na(phenotypes_nonan)] <- FALSE

phecodes <- data.table(phecodes = names(phenotypes_nonan)[-1])
phecodes <- addPhecodeInfo(phecodes)

## Processing auxiliary phewas data and genotype of interest

gene_cols <- colnames(genodata)[6:ncol(genodata)]

pb <- progress_bar$new(total = length(gene_cols), format = "[:bar] :percent :eta")

phe_results <- list()

covariates <- cbind(id.sex, genodata[, c("age")], genodata[, c("BMI")])
colnames(covariates)[3] <- "age"
colnames(covariates)[4] <- "BMI"
covariates$patient_id <- as.character(covariates$patient_id)

## Running PheWAS

for (gene in gene_cols) {
  pb$tick()
  gene <- gsub("-", ".", gene)
  gene <- gsub(":", ".", gene)
  print(gene)
  genotype <- genodata[, c("patient_id", gene)]

  joined_genos <- inner_join(phenotypes_nonan, genotype, by = "patient_id")
  
  covariates <- cbind(id.sex, genodata[, c("age")], genodata[, c("BMI")])
  covariates.filtered <- filter(covariates, patient_id %in% joined_genos$patient_id)
  
  cov_gen_data = inner_join(inner_join(covariates.filtered, genotype),
                            phenotypes_nonan)
  colnames(cov_gen_data)[3] <- "age"
  colnames(cov_gen_data)[4] <- "BMI"
  cov_gen_data$patient_id <- as.character(cov_gen_data$patient_id)
  
  results <- phewas(phenotypes=names(phenotypes_nonan)[-1], genotypes=c(gene), 
                    covariates=c("age", "gender", "BMI"), data=cov_gen_data, 
                    cores=4, significance.threshold = c("fdr"))
  
  
  phe_results[[gene]] <- results
}

## Saving results

saveRDS(phe_results, file = "outputs/phewas_results_age_gender_bmi_adjusted.rds")
