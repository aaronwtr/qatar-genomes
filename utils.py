import os
from tqdm import tqdm
import pandas as pd
from scipy.stats import fisher_exact
import numpy as np


def extract_patient_phenotypes(df):
    """
    This function takes in the Qatari genome data and creates a dictionary with an entry for each patient, where the key
    is the patient ID and the value is a list of the patient's phenotypes. Note that the phenotypes are currently kept
    in a single string, separated by commas.
    """
    patient_phenotypes = {}
    for index, row in df.iterrows():
        patient_id = row[os.getenv("PATIENT_ID")]
        phenotype = str(row[os.getenv("PHENOTYPE_ID")])
        if phenotype == 'nan' or patient_id == np.nan:
            continue
        phenotypes = row[os.getenv("PHENOTYPE_ID")].split(';')
        patient_phenotypes[int(patient_id)] = phenotypes
    return patient_phenotypes


def count_diabetic_patients(patient_dict):
    """
    This function takes in the patient dictionary and counts the number of patients with diabetes. Care must be taken
    because there is no notational convention of the diabetes phenotype in the data. The function returns a dictionary
    with the number of diabetic patients and the number of non-diabetic patients.
    """
    diabetic_dict = {}
    diabetes_entry = 0
    for patient_id, phenotypes in tqdm(patient_dict.items()):
        for phenotype in phenotypes:
            if 'diabetes' in phenotype.lower() and 'type 1' not in phenotype.lower():
                diabetes_entry += 1
        if diabetes_entry > 0:
            diabetes_entry = 0
            diabetic_dict[patient_id] = 'diabetic'
        else:
            diabetic_dict[patient_id] = 'non-diabetic'
    return diabetic_dict


def count_gene_hets_homs(gene_data, diabetes_status):
    """
    This function takes in the gene name and returns a dictionary with the number of heterozygous and homozygous
    patients for that gene.
    """
    gene = gene_data.columns[-1]
    gene_data = gene_data.dropna(how='all')
    gene_data = gene_data.replace("Hom", "Hom_mut")
    gene_data = gene_data.fillna("Hom_wt")
    print(gene_data)
    # TODO 1: Inspect gene_data for similarity with excel sheet. We should find 6 hom_mut diabetics
    patient_alleles = {}
    for index, row in gene_data.iterrows():
        patient_id = row[os.getenv("PATIENT_ID")]
        allele = row[gene]
        patient_alleles[patient_id] = allele
    patient_alleles_diabetes = {}
    for patient_id, allele in patient_alleles.items():
        patient_alleles_diabetes[patient_id] = [allele, diabetes_status[patient_id]]

    hom_mut_diab, hom_wt_diab, het_diab = 0, 0, 0
    hom_mut_non_diab, hom_wt_non_diab, het_non_diab = 0, 0, 0
    for patient_id, allele_diabetes in patient_alleles_diabetes.items():
        if allele_diabetes[0] == 'Hom_mut' and allele_diabetes[1] == 'diabetic':
            hom_mut_diab += 1
        elif allele_diabetes[0] == 'Hom_wt' and allele_diabetes[1] == 'diabetic':
            hom_wt_diab += 1
        elif allele_diabetes[0] == 'Het' and allele_diabetes[1] == 'diabetic':
            het_diab += 1
        elif allele_diabetes[0] == 'Hom_mut' and allele_diabetes[1] == 'non-diabetic':
            hom_mut_non_diab += 1
        elif allele_diabetes[0] == 'Hom_wt' and allele_diabetes[1] == 'non-diabetic':
            hom_wt_non_diab += 1
        elif allele_diabetes[0] == 'Het' and allele_diabetes[1] == 'non-diabetic':
            het_non_diab += 1
    return [hom_mut_diab, hom_wt_diab, het_diab, hom_mut_non_diab, hom_wt_non_diab, het_non_diab]


def make_cont_table(allele_counts):
    """
    This function takes in the allele counts and returns a contingency table.
    """
    cont_table = pd.DataFrame(
        [[allele_counts[0], allele_counts[1], allele_counts[2]],
         [allele_counts[3], allele_counts[4], allele_counts[5]]],
        index=['Diabetic', 'Non-diabetic'], columns=['Hom_mut', 'Hom_wt', 'Het'])
    cont_table_hom_mut = cont_table[['Hom_wt', 'Hom_mut']].transpose()
    cont_table_het = cont_table[['Hom_wt', 'Het']].transpose()
    return cont_table_hom_mut, cont_table_het


def fisher_exact_test(cont_table):
    """
    This function takes in a contingency table and returns the p-value of the Fisher's exact test.
    """
    if type(cont_table) == pd.DataFrame:
        cont_table = cont_table.to_numpy()
    _, pvalue = fisher_exact(cont_table)
    return pvalue
