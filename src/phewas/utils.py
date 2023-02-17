from tqdm import tqdm
import pickle as pkl
import pandas as pd
import numpy as np
import os


def open_icd10_mapping(path_to_map):
    """Open the ICD10 mapping file and return a dictionary of mappings."""
    with open(path_to_map, 'r') as f:
        mapping = {}
        for line in f:
            line = line.strip().split(',')
            icd10 = line[0]
            phenotype = line[1][0:].lower()
            mapping[phenotype] = icd10
    return mapping


def open_phecode_mapping(path_to_map):
    """Open the ICD10 mapping file and return a dictionary of phecode to ICD10 mappings."""
    with open(path_to_map, 'r') as f:
        mapping = {}
        for line in f:
            line = line.strip().split(',')
            icd10 = line[0]
            phenotype = line[3][0:].lower()
            mapping[phenotype] = icd10
    return mapping


def open_snomed_mapping(path_to_map, icd10_map):
    """Open the SNOMED mapping file and return a dictionary of mappings."""
    with open(path_to_map, 'r') as f:
        mapping = {}
        for line in f:
            line = line.strip().split('\t')
            phenotype = line[7].lower()
            try:
                icd10 = icd10_map[phenotype]
            except KeyError:
                continue
            mapping[phenotype] = icd10
    return mapping


def patient_icd10_map(patient_phenotype_dict, icd10_map, phecode_map, snomed_map):
    """Map the phenotypes to ICD10 and save in a dictionary where the key is the patient id and the values are the
    ICD10 codes for that patient."""
    if os.path.exists('../data/full_icd10_map.pkl'):
        with open('../data/full_icd10_map.pkl', 'rb') as f:
            merged = pkl.load(f)
        f.close()

        with open('../data/no_icd10_found.pkl', 'rb') as f:
            no_icd10_found = pkl.load(f)
        f.close()
        return merged, no_icd10_found
    else:
        patient_icd10 = {}
        no_icd10_found = {}
        with open('../data/manual_icd10_entries.csv', 'r') as f:
            manual_icd10 = {}
            for line in f:
                if line.startswith('\ufeffterm'):
                    continue
                line = line.strip().split(',')
                manual_icd10[line[0]] = line[1]
        for patient, phenotypes in tqdm(patient_phenotype_dict.items()):
            icd10_codes = []
            for phenotype in phenotypes:
                if phenotype in list(icd10_map.keys()):
                    if ' ' not in icd10_map[phenotype]:
                        icd10_codes.append(icd10_map[phenotype])
                elif phenotype in list(phecode_map.keys()):
                    if ' ' not in phecode_map[phenotype]:
                        icd10_codes.append(phecode_map[phenotype])
                elif phenotype in list(snomed_map.keys()):
                    if ' ' not in snomed_map[phenotype]:
                        icd10_codes.append(snomed_map[phenotype])
                else:
                    if patient not in no_icd10_found.keys():
                        no_icd10_found[patient] = [phenotype]
                    else:
                        no_icd10_found[patient].append(phenotype)
            patient_icd10[patient] = icd10_codes
        # TODO Fix this. Something is wrong with the manual icd10 entries.
        merged = {**manual_icd10, **patient_icd10}
        print(merged)

        with open('../data/full_icd10_map.pkl', 'wb') as f:
            pkl.dump(merged, f)
        f.close()

        with open('../data/no_icd10_found.pkl', 'wb') as f:
            pkl.dump(no_icd10_found, f)

    return merged, no_icd10_found


def calculate_unique_phenotypes(patient_icd10_dict):
    """Calculate the number of unique phenotypes in the patient_icd10 dictionary."""
    unique_phenotypes = set()
    for key, values in patient_icd10_dict.items():
        for value in values:
            unique_phenotypes.add(value)
    return len(unique_phenotypes)


def make_phewas_table(qatari_data, patient_phenotypes, gene):
    phewas_df = qatari_data.loc[:, ['Dummy ID for GEL', 'age', gene, 'gender']]
    phewas_df = phewas_df.rename(columns={'Dummy ID for GEL': 'patient_id'})
    phewas_df[gene] = phewas_df[gene].fillna(0)
    phewas_df[gene] = phewas_df[gene].replace('Hom', 1)
    phewas_df[gene] = phewas_df[gene].replace('Het', 2)
    unique_phenotypes = set()
    for key, values in patient_phenotypes.items():
        for value in values:
            unique_phenotypes.add(value)
    unique_phenotypes = list(unique_phenotypes)
    patient_ids = phewas_df['patient_id'].tolist()
    phewas_array = phewas_df['patient_id'].to_numpy().reshape(-1, 1)
    # constructing phenotype matrix for phewas
    for phenotype in tqdm(unique_phenotypes):
        phenotype_column = []
        for patient_id in patient_ids:
            if phenotype in patient_phenotypes[patient_id]:
                phenotype_column.append(1)
            else:
                phenotype_column.append(0)
        phenotype_column = np.array(phenotype_column).reshape(-1, 1)
        phewas_array = np.hstack((phewas_array, phenotype_column))

    phen_df = pd.DataFrame(phewas_array)
    phen_df = phen_df.drop(0, axis=1)
    phen_df.columns = unique_phenotypes
    phewas_df = pd.concat([phewas_df, phen_df], axis=1)
    print(phewas_df)

