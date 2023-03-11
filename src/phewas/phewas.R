library(PheWAS)
library(data.table)
library(tidyverse)
library(readxl)


setwd("~/Desktop/PhD/Research/QMUL/Research/Qatar Genomes Project/qatar-genomes")
#(data <- read_excel("C:/Enserio/OneDrive - Queen Mary, University of London/Projects/AWenteler/TRPV1_phewas_table.xlsx",sheet = "Sheet1"))

data <- fread("data/phewas_tables/TRPV1_phewas_table.csv")

## Processing phenotype data

# Create sex table
id.sex <- data[,c("patient_id","gender")]
id.sex$gender[which(id.sex$gender=="MALE")]<-"M"
id.sex$gender[which(id.sex$gender=="FEMALE")]<-"F"

data$patient_id <- as.character(data$patient_id)
data$gender <- as.factor(data$gender)
id.sex$patient_id <- as.character(id.sex$patient_id)
id.sex$gender <- as.factor(id.sex$gender)

# Read pheno test table
phenodata <- fread("data/phewas_tables/TRPV1_icd10_test.csv")

# Create sex table matching the test phenotype table (test subset)
id.sex.test <- id.sex[which(id.sex$patient_id %in% as.character(phenodata$id)),]
head(id.sex.test)

colnames(phenodata)
colnames(id.sex)

colnames(phenodata)[1] <- "patient_id"

# TODO Add full.population.ids argument
phenotypes <- createPhenotypes(phenodata, min.code.count = 1, 
                               id.sex = id.sex.test, 
                               vocabulary.map = PheWAS::phecode_map_icd10)

## Processing auxiliary phewas data and genotype of interest
genotype <- data[, c("TRPV1")]

# TODO Download the example phewas script to inspect how covariates should look like
covariates <- cbind(id.sex, data[, c("age")])




