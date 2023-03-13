library(PheWAS)
library(data.table)
library(tidyverse)
library(readxl)
library(dplyr)


setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")
#(data <- read_excel("C:/Enserio/OneDrive - Queen Mary, University of London/Projects/AWenteler/TRPV1_phewas_table.xlsx",sheet = "Sheet1"))

data <- fread("data/phewas_tables/TRPV1_phewas_table.csv")

## Processing phenotype data

id.sex <- data[,c("patient_id","gender")]
id.sex$gender[which(id.sex$gender=="MALE")]<-"M"
id.sex$gender[which(id.sex$gender=="FEMALE")]<-"F"

data$patient_id <- as.character(data$patient_id)
data$gender <- as.factor(data$gender)
id.sex$patient_id <- as.character(id.sex$patient_id)
id.sex$gender <- as.factor(id.sex$gender)
id.sex$ismale <- ifelse(id.sex$gender == 'M', TRUE, FALSE)
id.sex <- subset(id.sex, select = -gender)

phenodata <- fread("data/phewas_tables/TRPV1_icd10_test.csv")

phenodata <- setDT(phenodata)
phenodata<- phenodata[, patient_id := as.character(patient_id)]
test_ids <- inner_join(phenodata, genotype, by = "patient_id")
id.sex.test <- filter(genotype, patient_id %in% joined_genos$patient_id)
head(id.sex.test)

colnames(phenodata)
colnames(id.sex)

colnames(phenodata)[1] <- "patient_id"

# TODO Add full.population.ids argument
phenotypes <- createPhenotypes(phenodata, min.code.count = 1, 
                               id.sex = id.sex.test, 
                               vocabulary.map = PheWAS::phecode_map_icd10)

## Processing auxiliary phewas data and genotype of interest
genotype <- data[, c("patient_id", "TRPV1")]
joined_genos <- inner_join(phenotypes, genotype, by = "patient_id")
genotype.test <- filter(genotype, patient_id %in% joined_genos$patient_id)

covariates <- cbind(id.sex, data[, c("age")])
covariates.test <- filter(covariates, patient_id %in% joined_genos$patient_id)

## Running PheWAS

test_data = inner_join(inner_join(covariates, genotype),phenotypes)
results=phewas(phenotypes=names(phenotypes)[-1], genotypes=c("TRPV1"), 
               covariates=c("age", "ismale"), data=test_data, cores=4)

save(results, "~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes/phewas.csv", 
     compress=FALSE)
