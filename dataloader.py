import pandas as pd


class DataLoader:
    def __init__(self, path):
        self.path = path
        self.df = pd.read_excel(self.path)

    def get_qatari_data(self):
        return self.df
