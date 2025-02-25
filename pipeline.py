import pandas as pd
import requests
import re
from utils import extract_hmr_id, extract_mar_id, fetch_bigg_from_api

class MetaFluxProcessor:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def process(self):
        df = pd.read_csv(self.input_file, header=None, names=["Raw", "Value"])
        df["MAR_ID"] = df["Raw"].apply(extract_mar_id)
        df.to_csv(self.output_file, index=False)

class HMRToHumanGEMConverter:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def convert(self):
        df = pd.read_csv(self.input_file)
        df["Human_GEM_ID"] = df["HMR_ID"].apply(extract_hmr_id)
        df.to_csv(self.output_file, index=False)

class BiGGIDFetcher:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def fetch_bigg_ids(self):
        df = pd.read_csv(self.input_file)
        df["BiGG_ID"] = df["MAR_ID"].apply(fetch_bigg_from_api)
        df.to_csv(self.output_file, index=False)
        
class ModuleAssigner:
    def __init__(self, input_file, module_file, output_file):
        self.input_file = input_file
        self.module_file = module_file
        self.output_file = output_file

    def assign_modules(self):
        df = pd.read_csv(self.input_file)
        df_modules = pd.read_csv(self.module_file)
        df = df.merge(df_modules, how="left", left_on="BiGG_ID", right_on="Gene").drop(columns=["Gene"])
        df_assigned = df.dropna(subset=["Module"])
        df_unassigned = df[df["Module"].isna()].drop(columns=["Module"])
        df_assigned.to_csv(self.output_file, index=False)
        df_unassigned.to_csv("Unassigned_Genes.csv", index=False)