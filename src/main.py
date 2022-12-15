from utils import *
from dataloader import DataLoader
from plot import bar_plot
from dotenv import load_dotenv
import os
import pickle as pkl


def find_diabetic_patients():
    print("Loading data...")
    qatari_data = DataLoader(os.getenv("DATA")).get_qatari_data()
    print("Data loaded.")
    print("Parsing patient phenotypes...")
    patient_phenotypes = extract_patient_phenotypes(qatari_data)
    print("Patient phenotypes parsed.")
    print("Counting diabetic patients...")
    diabetic_patients = count_diabetic_patients(patient_phenotypes)
    print("Diabetic patients counted.")
    diabetic_bar_plot = bar_plot(diabetic_patients)
    return qatari_data, diabetic_patients


def het_hom_analysis(gene):
    qatari_data, diabetic_patients = find_diabetic_patients()
    gene_data = DataLoader(qatari_data).get_gene_data(gene)
    print("Counting alleles for different phenotypes...")
    counts = count_gene_hets_homs(gene_data, diabetic_patients)
    print("Alleles counted.")
    print("Making contingency table...")
    cont_table_hom_mut, cont_table_het = make_cont_table(counts)
    print(cont_table_hom_mut)
    print(cont_table_het)
    print("Contingency table made.")
    p_val_hom_mut = fisher_exact_test(cont_table_hom_mut)
    p_val_het = fisher_exact_test(cont_table_het)
    print("p-value for homozygote mutant alleles: ", p_val_hom_mut)
    print("p-value for heterozygote alleles: ", p_val_het)


if __name__ == "__main__":
    load_dotenv()
    het_hom_analysis("TRPV1")
