import pandas as pd


class DataLoader:
    def __init__(self, data):
        self.data = data
        if isinstance(data, str):
            self.datatype = data.split(".")[-1]
        elif isinstance(data, pd.DataFrame):
            self.datatype = "df"
        else:
            raise TypeError("Data type not supported. Choose either .xlsx, .pkl, or a pandas dataframe.")

    def get_qatari_data(self):
        if self.datatype == "xlsx":
            return pd.read_excel(self.data)
        elif self.datatype == "pkl":
            return pd.read_pickle(self.data)
        else:
            raise TypeError("Data type not supported. Choose either .xlsx or .pkl")

    def get_gene_data(self, gene):
        if self.datatype == "xlsx":
            df = pd.read_excel(self.data)
        elif self.datatype == "pkl":
            df = pd.read_pickle(self.data)
        elif self.datatype == "df":
            return self.data.iloc[:, 0:4].join(self.data[gene])
        else:
            raise TypeError("Data type not supported. Choose either .xlsx or .pkl")
        return df.iloc[:, 0:4].join(df[gene])
