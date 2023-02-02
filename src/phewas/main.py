from dotenv import load_dotenv
import os

from utils import *
from src.gwas.utils import splitted_patient_phenotypes
from dataloader import DataLoader


def phenotype_preprocessing():
    """
    This function should create a dictionary such that each key is a patient ID and a value is a list of homogenized
    phenotypes. Homogenization entails that irregularies in the phenotype strings are removed, such as spaces,
    capitalization, and punctuation.
    :return:
    """
    print("Loading data...")
    qatari_data = DataLoader(os.getenv("DATA")).get_qatari_data()
    print("Data loaded.")
    print("Parsing patient phenotypes...")
    patient_phenotypes = splitted_patient_phenotypes(qatari_data)
    print("Patient phenotypes parsed.")
    print("Performing Phenotype to ICD10 mapping...")
    icd10_mapping = open_icd10_mapping(os.getenv("ICD10_MAP"))
    phecode_mapping = open_phecode_mapping(os.getenv("ICD10_MAP"))
    patient_icd10, no_icd10_found = patient_icd10_map(patient_phenotypes, icd10_mapping, phecode_mapping)
    no_icd10 = []
    for key, values in no_icd10_found.items():
        for value in values:
            no_icd10.append(value)

    icd10 = []
    for key, values in patient_icd10.items():
        for value in values:
            icd10.append(value)

    print(len(icd10))
    print(len(no_icd10))


if __name__ == "__main__":
    load_dotenv()
    phenotype_preprocessing()
