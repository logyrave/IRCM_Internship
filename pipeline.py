import pandas as pd
import requests
import re
from utils import extract_hmr_id, extract_mar_id, batch_fetch_bigg, process_metadata_column, process_raw_Metaflux_output

class MetaFluxProcessor:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def process(self):
        df = process_raw_Metaflux_output(self.input_file)

        # 🔹 Assurer que MAR_ID est bien généré avant de passer à batch_fetch_bigg
        df["MAR_ID"] = df["Metadata"].apply(extract_mar_id)

        # 🔹 Supprimer les éventuelles lignes NaN introduites par l'opération
        df.dropna(subset=["Metadata"], inplace=True)

        df.to_csv(self.output_file, index=False)
        print(f"✅ Fichier {self.output_file} généré avec succès !")
        
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
        batch_fetch_bigg(self.input_file, self.output_file)
        df.to_csv(self.output_file, index=False)
        
class ModuleAssigner:
    def __init__(self, input_file, module_file, output_file):
        self.input_file = input_file
        self.module_file = module_file
        self.output_file = output_file
        
    def assign_modules(self):
        df = pd.read_csv(self.input_file)
        df = process_metadata_column(df)
        df_modules = pd.read_csv(self.module_file)
        df["BiGG_ID"] = df["BiGG_ID"].astype(str)
        df_modules["Gene"] = df_modules["Gene"].astype(str)
        df = df.merge(df_modules, how="left", left_on="BiGG_ID", right_on="Gene").drop(columns=["Gene"])
        df_assigned = df.dropna(subset=["Module"])
        df_unassigned = df[df["Module"].isna()].drop(columns=["Module"])
        df_assigned.to_csv(self.output_file, index=False)
        df_unassigned.to_csv("Unassigned_Genes.csv", index=False)
         
         
'''
add feature that will compare the output with scFEA 
+ add global statistics (e.g. number of genes assigned, number of genes unassigned, etc.)
'''