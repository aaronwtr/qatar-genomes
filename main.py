from utils import extract_patient_phenotypes, count_diabetic_patients, count_gene_hets_homs
from dataloader import DataLoader
from plot import bar_plot
from dotenv import load_dotenv
import os
import pickle as pkl


def entire_dataset_analysis():
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


def het_hom_analysis(gene):
    print("Loading data...")
    qatari_data = DataLoader(os.getenv("DATA")).get_gene_data(gene)
    print("Data loaded.")
    gene_data = count_gene_hets_homs(qatari_data)


if __name__ == "__main__":
    load_dotenv()
    het_hom_analysis("TRPV1")
