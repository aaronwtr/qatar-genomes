import pandas as pd
from dotenv import load_dotenv
import numpy as np
import nltk

from utils import *
from src.gwas.utils import splitted_patient_phenotypes, corpusify_phenotypes, tokenize_phewas_corpus, create_wordcloud, \
    map_tokens_to_phenotypes
from dataloader import DataLoader


def phenotype_preprocessing(wc=False, test=False):
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
    if test:
        transform_phewas_table(patient_icd10)
    mappings = {**icd10_mapping, **phecode_mapping, **snomed_mapping}
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
    if wc:
        corpus = corpusify_phenotypes(patient_icd10)
        tokenized_corpus = tokenize_phewas_corpus(corpus)
        phenotype_counts = map_tokens_to_phenotypes(tokenized_corpus, mappings)
        phenotype_counts = nltk.FreqDist(phenotype_counts)
        print("Patient phenotypes tokenized.")
        print("Creating wordcloud...")
        create_wordcloud(phenotype_counts, title="Wordcloud of mapped phenotypes")

        corpus = corpusify_phenotypes(no_icd10_found, mapped=False)
        tokenized_corpus = tokenize_phewas_corpus(corpus)
        phenotype_counts = dict(tokenized_corpus)
        create_wordcloud(phenotype_counts, title="Wordcloud of non-mapped phenotypes")

    return qatari_data


def phewas(patient_phenotypes, genes):
    """
    This function should create a dataframe that contains the covariates and phenotypes needed for a phewas analysis.
    This function also
    """
    print("Loading data...")
    qatari_data = DataLoader(os.getenv("DATA")).get_qatari_data()
    print("Data loaded.")
    print("Generating PheWAS table...")
    make_phewas_table(qatari_data, patient_phenotypes, genes)


if __name__ == "__main__":
    load_dotenv()
    files = os.listdir('../data')
    #if 'full_icd10_map.pkl' not in files:
    qatari_data = phenotype_preprocessing()
    genes = list(qatari_data.columns[4:])
    gene = "TRPV1"
    patient_phenotypes = pd.read_pickle('../data/full_icd10_map.pkl')
    phewas(patient_phenotypes, genes)

    # TODO: Write mapping to csv and look up the phenotypes to be inspected manually. Then prepare a csv with only these
    #  phenotypes and their ICD10 codes.

    # TODO: Adjust make_phewas_table such that it contains all the genes in our dataset, not just one.
