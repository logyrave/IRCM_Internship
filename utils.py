import re
import requests

def extract_mar_id(text):
    match = re.search(r"MAR0\d+", str(text))
    return match.group(0) if match else None

def extract_hmr_id(text):
    match = re.search(r"HMR_\d+", str(text))
    return match.group(0) if match else None

def fetch_bigg_from_api(mar_id):
    if not mar_id:
        return None
    url = f"https://metabolicatlas.org/api/v2/reactions/{mar_id}?model=HumanGem"
    try:
        response = requests.get(url, headers={"accept": "application/json"}, timeout=10)
        response.raise_for_status()
        data = response.json()
        bigg_entry = data.get("externalDbs", {}).get("BiGG", [])
        return bigg_entry[0]["id"] if bigg_entry else None
    except requests.exceptions.RequestException:
        return None
