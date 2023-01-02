import os
from tqdm import tqdm
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import numpy as np
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from nltk.probability import FreqDist
from wordcloud import WordCloud


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


def get_logreg_features(gene_data):
    gene_data = gene_data.dropna(how='all')
    gene_data = gene_data.replace("Hom", "Hom_mut")
    gene_data = gene_data.fillna("Hom_wt")
    gene_data['diabetes'] = gene_data['ICD-10 phenotype'].apply(lambda x: 1 if 'diabetes' in x.lower() and 'type 1'
                                                                               not in x.lower() else 0)
    gene_data['gender_bin'] = gene_data['gender'].apply(lambda x: 1 if x == 'FEMALE' else 0)
    gene_data['allele_hom_mut'] = gene_data[gene_data.columns[4]].apply(lambda x: 1 if x == 'Hom_mut' else 0)
    gene_data['allele_het'] = gene_data[gene_data.columns[4]].apply(lambda x: 1 if x == 'Het' else 0)
    features_hom_mut = gene_data[['gender_bin', 'age', 'allele_hom_mut', 'diabetes']]
    features_het = gene_data[['gender_bin', 'age', 'allele_het', 'diabetes']]
    return features_hom_mut, features_het


def fisher_exact_test(cont_table):
    """
    This function takes in a contingency table and returns the p-value of the Fisher's exact test.
    """
    if type(cont_table) == pd.DataFrame:
        cont_table = cont_table.to_numpy()
    _, pvalue = fisher_exact(cont_table)
    return round(pvalue, 3)


def corpusify_phenotypes(phenotypes):
    """
    This function takes in a dictionary of phenotypes and returns a corpus of the phenotypes.
    """
    corpus = ''
    phenotypes = list(phenotypes.values())
    for phenotype in phenotypes:
        corpus += phenotype[0] + ' '
    return corpus


def tokenize_corpus(corpus):
    """
    This function takes in a corpus and returns a list of tokens.
    """
    exclusion_list = ['low', 'check', 'foot', 'requested', 'current', 'persons', 'upper', 'examination',
                      'unspecified', 'complication', 'delivery', 'laboratory', 'cohort', 'consultation', 'prescription'
                      'classified', 'participate', 'single', 'test', 'explanation', 'investigation', 'suspected', 'use',
                      'right', 'left', 'history', 'essential', 'person', 'activity', 'prescription', 'telephone',
                      'without', 'specified', 'disease', 'screening', 'need', 'type', 'diseases', 'syndrome', 'vitamin',
                      'health', 'due', 'severe', 'results', 'findings', 'elsewhere', 'medical']
    stop_words = set(stopwords.words('english'))
    filtered_words = [word.lower() for word in word_tokenize(corpus) if word.lower() not in stop_words
                      and word.lower() not in exclusion_list]
    filtered_words = [word for word in filtered_words if word.isalpha()]
    return filtered_words


def create_wordcloud(tokens):
    """
    This function creates a wordcloud of the phenotypes.
    """
    fdist = FreqDist(tokens)
    wordcloud = WordCloud(width=800, height=800, background_color='white', min_font_size=10)\
        .generate_from_frequencies(fdist)
    plt.figure(figsize=(8, 8), facecolor=None)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.savefig('plots/wordcloud.png')

    return fdist


def get_icd10_matches(qatari_data, top_phenotype_tokens, n=100):
    qatari_phenotypes = extract_patient_phenotypes(qatari_data)
    qatari_phenotypes = [phenotype[0].split(',') for phenotype in qatari_phenotypes.values()]
    qatari_phenotypes_filt = list(set([phenotype for sublist in qatari_phenotypes for phenotype in sublist]))
    qatari_phenotypes_filt = [phenotype for phenotype in qatari_phenotypes_filt if phenotype != '']
    top_phenotype_tokens = top_phenotype_tokens.most_common(n)
    top_phenotype_tokens = dict(top_phenotype_tokens)
    top_tokens_list = list(top_phenotype_tokens.keys())
    top_icd10 = {}
    for phenotype in qatari_phenotypes_filt:
        for token in top_tokens_list:
            if token in phenotype:
                phe_count = top_phenotype_tokens[token]
                top_icd10[phenotype] = phe_count
    top_icd10 = dict(sorted(top_icd10.items(), key=lambda item: item[1], reverse=True))
    top_icd10_counted = {}
    for key in tqdm(top_icd10.keys()):
        top_icd10_counted[key] = sum(key in sublist for sublist in qatari_phenotypes)
    top_icd10_counted = dict(sorted(top_icd10_counted.items(), key=lambda item: item[1], reverse=True))
    with open('data/top_icd10_counted.csv', 'w') as f:
        for key in top_icd10_counted.keys():
            f.write("%s,%s\n" % (key, top_icd10_counted[key]))
    f.close()
    return top_icd10_counted
