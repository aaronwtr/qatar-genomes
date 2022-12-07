import numpy as np
import os


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
        phenotypes = row[os.getenv("PHENOTYPE_ID")].split(',')
        patient_phenotypes[int(patient_id)] = phenotypes
    return patient_phenotypes


def count_diabetic_patients(patient_dict):
    """
    This function takes in the patient dictionary and counts the number of patients with diabetes. Care must be taken
    because there is no notational convention of the diabetes phenotype in the data. The function returns a dictionary
    with the number of diabetic patients and the number of non-diabetic patients.
    """
    diabetic_dict = {}
    diabetic_patients = 0
    diabetes_entry = 0
    for patient_id, phenotypes in patient_dict.items():
        for phenotype in phenotypes:
            if 'diabetes' in phenotype.lower() and 'type 1' not in phenotype.lower():
                diabetes_entry += 1
        if diabetes_entry > 0:
            diabetic_patients += 1
            diabetes_entry = 0
    diabetic_dict['diabetic'] = diabetic_patients
    diabetic_dict['non-diabetic'] = len(patient_dict) - diabetic_patients
    return diabetic_dict
