import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import math
import re 
import os 
from scipy.stats import zscore, norm
import ace_tools_open as tools
import seaborn as sns

'''
This part as for purpose to compare the value of both algorithms by computing Z-score on cellType 
and analyse the biological relevance of both method
'''
Metaflux = pd.read_csv('crash_test.csv')
scFEA = pd.read_csv('Liver_Test_scFEA.csv')

'''
apply good cell type to the metaflux data
'''
Metaflux['Batch_Num'] = Metaflux['Batch'].str.extract(r'(\d+)').astype(int)
buffer = {}
counter = 1
for batch in sorted(Metaflux['Batch_Num'].unique()):  
    unique_types = Metaflux.loc[Metaflux['Batch_Num'] == batch, 'Cell_Type'].unique()  
    for cell_type in sorted(unique_types):  
        buffer[(batch, cell_type)] = counter
        counter += 1
Metaflux['Cell_Type'] = Metaflux.apply(lambda row: buffer[(row['Batch_Num'], row['Cell_Type'])], axis=1)
Metaflux.drop(columns=['Batch_Num'], inplace=True)

#compute Z-score by Cell_Type
Metaflux['Z_Score'] = Metaflux.groupby(['Batch', 'Cell_Type'])['Value'].transform(zscore)
print (f"\n---------\n Metaflux : \n {Metaflux}\n---------\n")

'''
attribute metadata to scFEA output
'''
metadata_file = "GSE115469_CellClusterType.txt"  
metadata = pd.read_csv(metadata_file, sep="\t") 
metadata = metadata[["CellName", "Cluster#", "CellType"]]
scFEA = scFEA.merge(metadata, left_on="Unnamed: 0", right_on="CellName", how="left")
scFEA.drop(columns=["CellName"], inplace=True)
#print (f"\n---------\n scFEA : \n {scFEA}\n---------\n")


'''
group scFEA by cluster and then compute Z-scorer (expected 20 lines)
'''

# apply Z score to scFEA datas 
df_clustered = scFEA.groupby("Cluster#").mean(numeric_only=True).reset_index()
df_clustered["CellType"] = scFEA.groupby("Cluster#")["CellType"].first().values
cols = ["Cluster#", "CellType"] + [col for col in df_clustered.columns if col not in ["Cluster#", "CellType"]]
df_clustered = df_clustered[cols]
numeric_cols = [col for col in df_clustered.columns if col not in ["Cluster#", "CellType"]]
df_clustered[numeric_cols] = df_clustered[numeric_cols].apply(zscore, axis=1)
scFEA = df_clustered
print (f"\n---------\n scFEA : \n {scFEA}\n---------\n")

'''
# Check Z-score normality 
plt.figure(figsize=(12, 6))
for index, row in scFEA.iloc[:, 2:].iterrows():
    sns.kdeplot(row, label=f"by cell {index}", alpha=0.5)
x = np.linspace(-3, 3, 100)
plt.plot(x, norm.pdf(x, 0, 1), 'k--', label="Normal standard")
plt.title("Z-Scores distribution")
plt.xlabel("Z-Score")
plt.ylabel("Density")
plt.legend()
plt.show()
'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test heatmap pls ignore go next section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
scfea_modules = [col for col in scFEA.columns if col not in ["Cluster#", "CellType"]]
df_metaflux_agg = Metaflux.groupby(["Cell_Type", "Module"])["Z_Score"].mean().reset_index()
df_metaflux_pivot = df_metaflux_agg.pivot(index="Cell_Type", columns="Module", values="Z_Score").fillna(0)
df_metaflux_pivot.index.name = "Cluster#"
df_metaflux_pivot.reset_index(inplace=True)
for module in scfea_modules:
    if module not in df_metaflux_pivot.columns:
        df_metaflux_pivot[module] = 0  # Ajouter les modules manquants avec des valeurs 0
df_metaflux_pivot = df_metaflux_pivot[["Cluster#"] + scfea_modules]
df_metaflux_pivot = df_metaflux_pivot.set_index("Cluster#").T
df_scfea_pivot = scFEA.set_index("Cluster#").drop(columns=["CellType"]).T
fig, axes = plt.subplots(1, 2, figsize=(14, 12), sharey=True)
sns.heatmap(df_scfea_pivot, cmap="coolwarm", center=0, ax=axes[0])
axes[0].set_title("scFEA - modules activity")
axes[0].set_xlabel("Cluster#")
axes[0].set_ylabel("Modules")
sns.heatmap(df_metaflux_pivot, cmap="coolwarm", center=0, ax=axes[1])
axes[1].set_title("METAFlux - modules activity (module completion with 0)")
axes[1].set_xlabel("Cluster#")
axes[1].set_ylabel("Modules")
plt.tight_layout()
plt.show()
'''
#%%%%%%%%%%%%%%%%%%%%%%%%%% scFEA vs METAFlux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scfea_modules = [col for col in scFEA.columns if col not in ["Cluster#", "CellType"]]
df_metaflux_agg = Metaflux.groupby(["Cell_Type", "Module"])["Z_Score"].mean().reset_index()
df_metaflux_pivot = df_metaflux_agg.pivot(index="Cell_Type", columns="Module", values="Z_Score").fillna(0)


common_modules = scFEA.columns.intersection(df_metaflux_pivot.columns)
df_scfea_filtered = scFEA.set_index("Cluster#").drop(columns=["CellType"])[common_modules].T
df_metaflux_filtered = df_metaflux_pivot[common_modules].T
celltypes = scFEA.set_index("Cluster#")["CellType"]  
df_scfea_filtered.columns = celltypes[df_scfea_filtered.columns].values
df_metaflux_filtered.columns = celltypes[df_metaflux_filtered.columns].values
rows = []
for module in common_modules:
    rows.append(df_scfea_filtered.loc[module]) 
    rows.append(df_metaflux_filtered.loc[module])  

# heatmap Z_score scFEA vs Z_score METAFlux (1 line by method alternatively)
df_combined = pd.DataFrame(rows)
df_combined.index = [f"{mod}_S" if i % 2 == 0 else f"{mod}_M" for i, mod in enumerate(common_modules.repeat(2))]
plt.figure(figsize=(14, 12))
cmap = sns.diverging_palette(10, 150, s=80, l=50, as_cmap=True)  # Rouge → Vert
sns.heatmap(df_combined, cmap=cmap, center=0, linewidths=0.5, linecolor="gray")
plt.title("scFEA vs METAFlux - Modules activity comparison")
plt.xlabel("CellType")
plt.ylabel("Modules (methods alternated)")
plt.xticks(rotation=90)
plt.show()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scFEA - METAFlux for direct comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
Heatmap scFEA (Z_score by cell type) - METAFlux (Z_score by cell type)
'''

# melt scFEA in long format to align with Metalfux
df_scfea_long = scFEA.melt(id_vars=['Cluster#', 'CellType'], 
                              var_name='Module', value_name='scFEA_Value')

# just in case 
df_metaflux_filtered = Metaflux[['Value', 'Module', 'Cell_Type', 'Z_Score']]
df_metaflux_filtered = df_metaflux_filtered.rename(columns={'Cell_Type': 'Cluster#', 'Value': 'METAFlux_Value'})


# type insurance
df_metaflux_filtered['Cluster#'] = df_metaflux_filtered['Cluster#'].astype(str)
df_scfea_long['Cluster#'] = df_scfea_long['Cluster#'].astype(str)

# merged on cluster x module
df_merged = pd.merge(df_scfea_long, df_metaflux_filtered, on=['Cluster#', 'Module'], how='left')

# pivot matrix
df_pivot = df_merged.pivot_table(index='Module', columns='CellType', values='METAFlux_Value')


# add column for =scFEA - METAFlux
df_merged['Diff_Value'] = df_merged['scFEA_Value'] - df_merged['METAFlux_Value']

# pivot matrix 
df_pivot_diff = df_merged.pivot_table(index='Module', columns='CellType', values='Diff_Value', fill_value=0)

# heatmap plot 
plt.figure(figsize=(12, 8))
sns.heatmap(df_pivot_diff, cmap="coolwarm", center=0, linewidths=0.5, linecolor="gray")

plt.title("Z_score scFEA - Z_score METAFlux ")
plt.xlabel("Cell_Type")
plt.ylabel("Module")
plt.xticks(rotation=45, ha="right")
plt.show()













