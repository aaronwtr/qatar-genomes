from utils import *
from dataloader import DataLoader
from plot import bar_plot
from dotenv import load_dotenv
import os


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


def fisher_allele_analysis(gene):
    """
    This function analyses the association between homozygous alleles and a particular phenotype, specifically diabetes.
    It transforms our data into a 2x2 contingency table to calculate the p-values to quantify the association of an
    allele with the phenotype.
    :param gene: Gene name you want to test for phenotypic association.
    :return: P-value for the allele and the phenotype association.
    """
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


def logreg_allele_analysis(gene):
    """
    This function analyses the association between alleles and a particular phenotype by using a simple logistic
    regression model.
    :param gene: Gene name you want to test for phenotypic association
    :return:
    """
    print("Getting feature data for logistic regression model...")
    qatari_data, diabetic_patients = find_diabetic_patients()
    gene_data = DataLoader(qatari_data).get_gene_data(gene)
    hom_mut_features, het_features = get_logreg_features(gene_data)
    hom_mut_model = smf.logit("diabetes ~ allele_hom_mut + age + gender_bin", hom_mut_features).fit()
    het_model = smf.logit("diabetes ~ allele_het + age + gender_bin", het_features).fit()
    print(f"Homologous mutant model: {hom_mut_model.summary()}")
    print(f"Heterozygous model: {het_model.summary()}")


if __name__ == "__main__":
    load_dotenv()
    logreg_allele_analysis("TRPV1")
