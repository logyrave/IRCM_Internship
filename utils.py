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
    """asynchrone request to fetch BiGG ID associated with MAR_ID."""
    
    url = API_URL.format(mar_id)
    try:
        async with session.get(url, headers={"accept": "application/json"}, ssl=False, timeout=10) as response:
            if response.status == 200:
                data = await response.json()
                print(f"API response for {mar_id}: {json.dumps(data, indent=2)}")  # DEBUG API

                # BiGG existence check :
                bigg_entry = data.get("externalDbs", {}).get("BiGG", [])
                return mar_id, bigg_entry[0]["id"] if bigg_entry else None
            else:
                print(f"⚠️ HTTP protocole error {response.status} for {mar_id}")
                return mar_id, None
    except Exception as e:
        print(f"❌ Fetch Error for {mar_id}: {e}")
        return mar_id, None

async def fetch_all_bigg(mar_ids):
    """parrarelise request by BiGG ID processed in batch."""
    async with aiohttp.ClientSession() as session:
        tasks = [fetch_bigg(session, mar_id) for mar_id in mar_ids]
        results = await asyncio.gather(*tasks)

        # trooble shooting
        print("\n Summary of API request result :")
        for mar, bigg in results[:10]:  
            print(f"MAR_ID: {mar} -> BiGG_ID: {bigg}")

        return dict(results)

def extract_hmr_to_mar(text):
    """Deal with both cases :
    1. ID MAR0xxxx → returned.
    2. ID HMR_xxxx → converted to MAR0xxxx.
    """
    text = str(text).strip()

    # Case where there is already a MAR0xxxx
    match_mar = re.search(r"(MAR0\d+)", text)
    if match_mar:
        return match_mar.group(1)

    # Converting HMR_xxxx to MAR0xxxx case
    match_hmr = re.search(r"HMR_(\d+)", text)
    if match_hmr:
        return f"MAR0{match_hmr.group(1)}"  

    # If no HMR / MAR ID identified then fuck off
    return None

def batch_fetch_bigg(input_file, output_file):
    """batch processing + ensure conversion ID."""
    df = pd.read_csv(input_file, low_memory=False)

    # check for header
    if "Metadata" in df.columns:
        print("⚠️ Header detected, first line passed.")
        df = df.iloc[1:].reset_index(drop=True)

    df["MAR_ID"] = df["Metadata"].apply(extract_hmr_to_mar)
    
    mar_ids = df["MAR_ID"].dropna().unique().tolist()
    bigg_mapping = asyncio.run(fetch_all_bigg(mar_ids))

    df["BiGG_ID"] = df["MAR_ID"].map(bigg_mapping)
    df.to_csv(output_file, index=False)
    print(df[["Metadata", "MAR_ID"]].drop_duplicates().head())  # check conversion
    print("BiGG ID fetching done")
    
    
def process_metadata_column(df):
    """First colonne separation in 1 colonne for each info"""
    df[["Batch", "Cell_Type"]] = df["Metadata"].str.extract(r"(batch_\d+)\.celltype (\d+)")
    df.drop(columns=["Metadata"], inplace=True)  
    return df

def process_raw_Metaflux_output(input_file):
    df = pd.read_csv(input_file, header=None, low_memory=False)

    # type / dim check
    if df.shape[1] < 2:
        raise ValueError("1 columns detected, at least 2 expected. (WTF?)")

    df.columns = ["Metadata"] + [f"V{i}" for i in range(1, df.shape[1])]

    df.iloc[:, 1:] = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")

    df["Metadata"] = df["Metadata"].apply(lambda x: re.sub(r"HMR_(\d+)", r"MAR0\1", str(x)))

    df["Value"] = df.iloc[:, 1:].mean(axis=1)

    # get the NaN out of the way
    df.dropna(subset=["Metadata"], inplace=True)

    df_final = df[["Metadata", "Value"]].copy()
    print(f"\nTransformation sucessed ! Final format :\n{df_final.head()}\n")
    return df_final

