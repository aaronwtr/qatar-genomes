from dotenv import load_dotenv
import os
import numpy as np

from utils import *
from src.gwas.utils import splitted_patient_phenotypes
from dataloader import DataLoader


def phenotype_preprocessing():
    """
    This function should create a dictionary such that each key is a patient ID and a value is a list of homogenized
    phenotypes. Homogenization entails that irregularities in the phenotype strings are removed, such as spaces,
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
    snomed_mapping = open_snomed_mapping(os.getenv("SNOMED_MAP"), icd10_mapping)
    patient_icd10, no_icd10_found = patient_icd10_map(patient_phenotypes, icd10_mapping, phecode_mapping,
                                                      snomed_mapping)
    print(icd10_mapping)
    unique_phens_found = calculate_unique_phenotypes(patient_icd10)
    unique_phens_not_found = calculate_unique_phenotypes(no_icd10_found)
    print(f"Mapped {unique_phens_found} unique phenotypes out of {unique_phens_not_found + unique_phens_found} ("
          f"{np.round(((unique_phens_found / (unique_phens_not_found + unique_phens_found)) * 100), 2)}%) total unique "
          f"phenotypes in the Qatari dataset.")
    no_icd10 = []
    for key, values in no_icd10_found.items():
        for value in values:
            no_icd10.append(value)

    icd10 = []
    for key, values in patient_icd10.items():
        for value in values:
            icd10.append(value)

    print(f"Mapped {len(icd10)} patient phenotypes out of {len(icd10) + len(no_icd10)} "
          f"({np.round((len(icd10) / (len(icd10) + len(no_icd10)) * 100), 2)}%) patient phenotypes in the Qatari "
          f"dataset.")


if __name__ == "__main__":
    load_dotenv()
    phenotype_preprocessing()
