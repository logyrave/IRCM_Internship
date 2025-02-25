import requests
import aiohttp
import asyncio
import pandas as pd
import re
import json 

def extract_mar_id(text):
    match = re.search(r"MAR0\d+", str(text))
    return match.group(0) if match else None

def extract_hmr_id(text):
    match = re.search(r"HMR_\d+", str(text))
    return match.group(0) if match else None


'''
Get BiGG gene ID from HMR reaction ID following the process: 
1. Extract MAR ID from HMR ID via simple regex
2. Request MetabolicAtlas API to get BiGG ID from MAR ID

async def used in order to parallelize requests because the API can take some time to respond
'''
API_URL = "https://metabolicatlas.org/api/v2/reactions/{}?model=HumanGem"

async def fetch_bigg(session, mar_id):
    """Effectue une requête asynchrone pour récupérer l'ID BiGG d'un MAR_ID."""
    
    url = API_URL.format(mar_id)
    try:
        async with session.get(url, headers={"accept": "application/json"}, ssl=False, timeout=10) as response:
            if response.status == 200:
                data = await response.json()
                print(f"✅ Réponse API pour {mar_id}: {json.dumps(data, indent=2)}")  # DEBUG API

                # Vérifier si BiGG existe dans la réponse JSON
                bigg_entry = data.get("externalDbs", {}).get("BiGG", [])
                return mar_id, bigg_entry[0]["id"] if bigg_entry else None
            else:
                print(f"⚠️ Erreur HTTP {response.status} pour {mar_id}")
                return mar_id, None
    except Exception as e:
        print(f"❌ Erreur lors de la récupération de {mar_id}: {e}")
        return mar_id, None

async def fetch_all_bigg(mar_ids):
    """Gère plusieurs requêtes asynchrones pour récupérer les IDs BiGG en batch."""
    async with aiohttp.ClientSession() as session:
        tasks = [fetch_bigg(session, mar_id) for mar_id in mar_ids]
        results = await asyncio.gather(*tasks)

        # Afficher les résultats pour voir s'ils sont bien extraits
        print("\n✅ Résumé des résultats BiGG récupérés :")
        for mar, bigg in results[:10]:  # Affiche les 10 premiers résultats
            print(f"MAR_ID: {mar} -> BiGG_ID: {bigg}")

        return dict(results)

def extract_hmr_to_mar(text):
    """Gère les cas où Metadata contient :
    1. Un ID MAR0xxxx → on le retourne directement.
    2. Un ID HMR_xxxx → on le convertit en MAR0xxxx.
    """
    text = str(text).strip()

    # Cas où Metadata contient déjà MAR0xxxx
    match_mar = re.search(r"(MAR0\d+)", text)
    if match_mar:
        return match_mar.group(1)

    # Cas où Metadata contient HMR_xxxx → conversion en MAR0xxxx
    match_hmr = re.search(r"HMR_(\d+)", text)
    if match_hmr:
        return f"MAR0{match_hmr.group(1)}"  # Conversion en MAR0xxxx

    # Si aucun des deux cas ne correspond, retourne None
    return None

def batch_fetch_bigg(input_file, output_file):
    """Lecture des données, conversion HMR -> MAR, et récupération des BiGG_ID en batch."""
    df = pd.read_csv(input_file, header=None, names=["Metadata", "Value"])
    df["MAR_ID"] = df["Metadata"].apply(extract_hmr_to_mar)
    
    mar_ids = df["MAR_ID"].dropna().unique().tolist()
    bigg_mapping = asyncio.run(fetch_all_bigg(mar_ids))
    
    df["BiGG_ID"] = df["MAR_ID"].map(bigg_mapping)
    df.to_csv(output_file, index=False)
    print(df[["Metadata", "MAR_ID"]].drop_duplicates().head())  # Vérifie conversion
    print("Récupération des BiGG_ID terminée.")
    
    
def process_metadata_column(df):
    """Sépare la colonne 'Metadata' en 'Batch' et 'Cell_Type' et supprime MAR_ID."""
    df[["Batch", "Cell_Type"]] = df["Metadata"].str.extract(r"(batch_\d+)\.celltype (\d+)")
    df.drop(columns=["Metadata"], inplace=True)  # Supprime la colonne originale
    return df