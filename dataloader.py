import pandas as pd


class DataLoader:
    def __init__(self, path):
        self.path = path
        self.datatype = path.split(".")[-1]

    def get_qatari_data(self):
        if self.datatype == "xlsx":
            return pd.read_excel(self.path)
        elif self.datatype == "pkl":
            return pd.read_pickle(self.path)
        else:
            raise TypeError("Data type not supported. Choose either .xlsx or .pkl")

    def get_gene_data(self, gene):
        if self.datatype == "xlsx":
            df = pd.read_excel(self.path)
        elif self.datatype == "pkl":
            df = pd.read_pickle(self.path)
        else:
            raise TypeError("Data type not supported. Choose either .xlsx or .pkl")
        return df.iloc[:, 0:4].join(df[gene])
