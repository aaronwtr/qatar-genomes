library(PheWAS)
library(data.table)
library(tidyverse)
library(readxl)
library(dplyr)
library(progress)


setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")

genodata <- fread("data/phewas_tables/gene_data.csv")
phenodata <- fread("data/phewas_tables/phenotype_data.csv")

## Processing phenotype data

id.sex <- genodata[,c("patient_id","gender")]
id.sex$gender[which(id.sex$gender=="MALE")]<-"M"
id.sex$gender[which(id.sex$gender=="FEMALE")]<-"F"

genodata$patient_id <- as.character(genodata$patient_id)
genodata$gender <- as.factor(genodata$gender)
id.sex$patient_id <- as.character(id.sex$patient_id)
id.sex$gender <- as.factor(id.sex$gender)

phenodata <- setDT(phenodata)
phenodata<- phenodata[, patient_id := as.character(patient_id)]
id.sex.test <- filter(id.sex, patient_id %in% phenodata$patient_id)
head(id.sex.test)

# TODO Add full.population.ids argument
phenotypes <- createPhenotypes(phenodata, min.code.count = 1, 
                               id.sex = id.sex.test, 
                               vocabulary.map = PheWAS::phecode_map_icd10,
                               rollup.map=PheWAS::phecode_rollup_map)
phenotypes_nonan<- phenotypes
phenotypes_nonan[is.na(phenotypes_nonan)]<-FALSE

## Processing auxiliary phewas data and genotype of interest

gene_cols <- colnames(genodata)[5:ncol(genodata)]
#genodata <- lapply(genodata,as.factor)
#genodata$patient_id <- as.numeric(as.character(genodata$patient_id))
#phenotypes$patient_id <- as.numeric(as.character(phenotypes$patient_id))
#genodata$`ICD-10 phenotype` <- as.character(genodata$`ICD-10 phenotype`)
#genodata$age <- as.numeric(genodata$age)
genodata <- data.frame(genodata)


pb <- progress_bar$new(total = length(gene_cols), format = "[:bar] :percent :eta")

phe_results <- list()

covariates <- cbind(id.sex, genodata[, c("age")])
colnames(covariates)[3] <- "age"
covariates$patient_id <- as.character(covariates$patient_id)

# all_data = inner_join(inner_join(covariates, genotype),phenotypes)


## Running PheWAS

for (gene in gene_cols) {
  pb$tick()
  gene <- gsub("-", ".", gene)
  genotype <- genodata[, c("patient_id", gene)]

  joined_genos <- inner_join(phenotypes_nonan, genotype, by = "patient_id")
  
  covariates <- cbind(id.sex, genodata[, c("age")])
  #covariates$patient_id <- as.numeric(as.character(covariates$patient_id))
  covariates.filtered <- filter(covariates, patient_id %in% joined_genos$patient_id)
  
  cov_gen_data = inner_join(inner_join(covariates.filtered, genotype),
                            phenotypes_nonan)
  colnames(cov_gen_data)[3] <- "age"
  cov_gen_data$patient_id <- as.character(cov_gen_data$patient_id)
  
  
  # cov_gen_data$gender <- as.factor(id.sex$gender)
  
  results <- phewas(phenotypes=names(phenotypes_nonan)[-1], genotypes=c(gene), 
                    covariates=c("age", "gender"), data=cov_gen_data, cores=4, 
                    significance.threshold = c("fdr"))
  
  
  phe_results[[gene]] <- results
}

## Saving results

output_path <- "outputs/phewas_full_not_unique.csv"

if (!file.exists(output_path)) {
  write.csv(results, file = output_path)
}