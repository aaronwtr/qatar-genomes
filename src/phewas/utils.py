from tqdm import tqdm
import pickle as pkl
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


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
    """
    Map the phenotypes to ICD10 and save in a dictionary where the key is the patient id and the values are the
    ICD10 codes for that patient.
    """
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
                elif phenotype in list(manual_icd10.keys()):
                    icd10_codes.append(manual_icd10[phenotype])
                else:
                    if patient not in no_icd10_found.keys():
                        no_icd10_found[patient] = [phenotype]
                    else:
                        no_icd10_found[patient].append(phenotype)
            patient_icd10[patient] = icd10_codes

        with open('../data/full_icd10_map.pkl', 'wb') as f:
            pkl.dump(patient_icd10, f)
        f.close()

        with open('../data/no_icd10_found.pkl', 'wb') as f:
            pkl.dump(no_icd10_found, f)

    return patient_icd10, no_icd10_found


def calculate_unique_phenotypes(patient_icd10_dict):
    """Calculate the number of unique phenotypes in the patient_icd10 dictionary."""
    unique_phenotypes = set()
    for key, values in patient_icd10_dict.items():
        for value in values:
            unique_phenotypes.add(value)
    return len(unique_phenotypes)


def make_phewas_table(qatari_data, patient_phenotypes):
    bmi_df = pd.read_excel('../data/phewas_input/BMI.xlsx')
    bmi_df = bmi_df.rename(columns={'Dummy ID for GEL': 'patient_id'})
    bmi_df = bmi_df[['patient_id', 'BodyFat - BMI']]
    phewas_df = qatari_data.rename(columns={'Dummy ID for GEL': 'patient_id'})
    gene_data = phewas_df.copy()
    gene_data.replace({np.nan: 0, 'Hom': 1, 'Het': 2}, inplace=True)
    gene_data_left = gene_data[['patient_id', 'ICD-10 phenotype', 'gender', 'age']]
    gene_data_right = gene_data.iloc[:, [0] + list(range(4, len(gene_data.columns)))]
    merged_df = pd.merge(gene_data_left, bmi_df, on='patient_id', how='left')
    gene_data = pd.merge(merged_df, gene_data_right, on='patient_id', how='left')

    phenotype_data = np.array([[key, value] for key in patient_phenotypes for value in patient_phenotypes[key]])
    code_col = np.full((phenotype_data.shape[0], 1), 'ICD10')
    reshape_phendat = np.empty((phenotype_data.shape[0], phenotype_data.shape[1] + 1), dtype=phenotype_data.dtype)
    reshape_phendat[:, 0] = phenotype_data[:, 0]
    reshape_phendat[:, 1] = code_col[:, 0]
    reshape_phendat[:, 2:] = phenotype_data[:, 1:]
    phenotype_data = reshape_phendat

    count_col = np.full((phenotype_data.shape[0], 1), 5)
    count_col[0] = 6
    phenotype_data = np.hstack((phenotype_data, count_col))

    phenotype_data = pd.DataFrame(phenotype_data)
    phenotype_data.columns = ['patient_id', 'vocabulary_id', 'code', 'count']

    phenotype_data.to_csv('../data/phewas_tables/phenotype_data.csv', index=False)
    gene_data.to_csv('../data/phewas_tables/gene_data_bmi.csv', index=False)

    return gene_data, phenotype_data


def transform_phewas_table(patient_icd10):
    manual_icd10 = pd.read_excel('../data/phewas_tables/TRPV1_phewas_table.xlsx', sheet_name=1)
    icd10_codes = list(manual_icd10.columns)[4:]
    patient_icd10 = {key: value for key, value in patient_icd10.items() if any(x in value for x in icd10_codes)}
    for key, value in patient_icd10.items():
        patient_icd10[key] = [x for x in value if x in icd10_codes]

    patient_ids = []
    icd10_codes = []
    vocab = []
    counts = []
    cnt = 0
    for key, value in patient_icd10.items():
        for code in value:
            patient_ids.append(key)
            icd10_codes.append(code)
            vocab.append('ICD10')
            if cnt == 0:
                counts.append(4)
            else:
                counts.append(3)
            cnt+=1
    patient_ids = np.array(patient_ids).reshape(-1, 1)
    icd10_codes = np.array(icd10_codes).reshape(-1, 1)
    vocab = np.array(vocab).reshape(-1, 1)
    counts = np.array(counts).reshape(-1, 1)
    icd10_array = np.hstack((patient_ids, vocab, icd10_codes, counts))
    icd10_array = pd.DataFrame(icd10_array)
    icd10_array.columns = ['id', 'vocabulary_id', 'code', 'count']
    icd10_array.to_csv('../data/phewas_tables/TRPV1_icd10_test.csv', index=False)


def count_pheno_groups(phenotypes):
    phenos = list(phenotypes['group'])
    pheno_count = {}
    for pheno in phenos:
        if pheno in pheno_count.keys():
            pheno_count[pheno] += 1
        else:
            pheno_count[pheno] = 1
    sorted_pheno_count = {k: v for k, v in sorted(pheno_count.items(), key=lambda item: item[1], reverse=True)}
    plt.figure(figsize=(10, 6))
    plt.bar(sorted_pheno_count.keys(), sorted_pheno_count.values())
    plt.xticks(rotation=45, fontsize=8)
    plt.tight_layout()
    plt.show()
    return pheno_count
