import requests
import aiohttp
import asyncio
import pandas as pd
import re

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
        async with session.get(url, headers={"accept": "application/json"}, timeout=10) as response:
            response.raise_for_status()
            data = await response.json()
            bigg_entry = data.get("externalDbs", {}).get("BiGG", [])
            return mar_id, bigg_entry[0]["id"] if bigg_entry else None
    except Exception as e:
        print(f"Erreur lors de la récupération de {mar_id}: {e}")
        return mar_id, None

async def fetch_all_bigg(mar_ids):
    """Gère plusieurs requêtes asynchrones pour récupérer les IDs BiGG en batch."""
    async with aiohttp.ClientSession() as session:
        tasks = [fetch_bigg(session, mar_id) for mar_id in mar_ids]
        results = await asyncio.gather(*tasks)
        return dict(results)

def extract_hmr_to_mar(text):
    """Convertit un identifiant HMR_xxxx en MAR0xxxx."""
    match = re.search(r"HMR_(\d+)", str(text))
    return f"MAR0{match.group(1)}" if match else None

def batch_fetch_bigg(input_file, output_file):
    """Lecture des données, conversion HMR -> MAR, et récupération des BiGG_ID en batch."""
    df = pd.read_csv(input_file)
    df["MAR_ID"] = df["Raw"].apply(extract_hmr_to_mar)
    
    mar_ids = df["MAR_ID"].dropna().unique().tolist()
    bigg_mapping = asyncio.run(fetch_all_bigg(mar_ids))
    
    df["BiGG_ID"] = df["MAR_ID"].map(bigg_mapping)
    df.to_csv(output_file, index=False)
    print("Récupération des BiGG_ID terminée.")