from tqdm import tqdm


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
    patient_icd10 = {}
    no_icd10_found = {}
    for patient, phenotypes in tqdm(patient_phenotype_dict.items()):
        icd10_codes = []
        for phenotype in phenotypes:
            if phenotype in list(icd10_map.keys()):
                icd10_codes.append(icd10_map[phenotype])
            elif phenotype in list(phecode_map.keys()):
                icd10_codes.append(phecode_map[phenotype])
            # elif phenotype in list(snomed_map.keys()):
            #     icd10_codes.append(snomed_map[phenotype])
            else:
                if patient not in no_icd10_found.keys():
                    no_icd10_found[patient] = [phenotype]
                else:
                    no_icd10_found[patient].append(phenotype)
        patient_icd10[patient] = icd10_codes
    return patient_icd10, no_icd10_found
