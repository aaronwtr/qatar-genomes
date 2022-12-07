from utils import extract_patient_phenotypes, count_diabetic_patients
from dataloader import DataLoader
from plot import bar_plot
from dotenv import load_dotenv
import os


def entire_dataset_analysis():
    qatari_data = DataLoader(os.getenv("DATA")).get_qatari_data()
    patient_phenotypes = extract_patient_phenotypes(qatari_data)
    diabetic_patients = count_diabetic_patients(patient_phenotypes)
    diabetic_bar_plot = bar_plot(diabetic_patients)


if __name__ == "__main__":
    load_dotenv()
    entire_dataset_analysis()
