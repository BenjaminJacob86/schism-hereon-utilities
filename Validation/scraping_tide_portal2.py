import requests
from bs4 import BeautifulSoup
import csv
import os
import time
import json

# Step 1: Fetch the main page (to allow any cookies/reference)
main_url = "https://www.kuestendaten.de/DE/Services/Messreihen_Dateien_Download/Download_Zeitreihen_node.html"
session = requests.Session()
resp = session.get(main_url)
resp.raise_for_status()

# Step 2: Inspect for XHR endpoint – assume we found it at `data_endpoint`
# (you’ll need to discover the actual endpoint from dev-tools)
data_endpoint = "https://www.kuestendaten.de/DE/…/stationlist.json"  # placeholder

resp2 = session.get(data_endpoint)
resp2.raise_for_status()
data = resp2.json()

# Example of what `data` might contain:
# [
#   {"stationName": "Abbenfleth Sperrwerk", "authority": "...", "linkId": "ABBENFLETH_SPERRWERK", "folder": "wiski", "parameter": ["Salzgehalt", "Wassertemperatur"]},
#   {"stationName": "D1 - Hanskalbsand, Oberfläche", "authority": "...", "linkId": "LZ_hanskalbsand_oberflaeche", "folder": "", "parameter": ["Salzgehalt"]},
#   ...
# ]

results = []
base = "https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download"

for rec in data:
    station = rec.get("stationName")
    authority = rec.get("authority")
    linkId = rec.get("linkId")
    folder = rec.get("folder", "")
    parameters = rec.get("parameter", [])
    for param in parameters:
        if folder:
            url = f"{base}/{folder}/{linkId}!{param}.zip"
        else:
            url = f"{base}/{linkId}!{param}.zip"
        results.append({"station": station, "authority": authority, "url": url, "parameter": param})

# Step 3: Save metadata
with open("stations_links.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["station","authority","parameter","url"])
    writer.writeheader()
    writer.writerows(results)

# Step 4: Download (optional)
for rec in results:
    url = rec["url"]
    fname = os.path.basename(url)
    print("Downloading", fname)
    try:
        d = session.get(url, stream=True)
        if d.status_code == 200:
            with open(fname, "wb") as f:
                for chunk in d.iter_content(chunk_size=8192):
                    f.write(chunk)
            print("  -> OK")
        else:
            print("  -> Failed status", d.status_code)
    except Exception as e:
        print("  -> Exception", e)
    time.sleep(0.5)

#-----------

## JSON based approach
import requests
import json
import csv
import os
import time

session = requests.Session()
base = "https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download"

# 1. Fetch JSONs
url_wiski = base + "/wiski/WISKI_datenDownload.json"
url_all = base + "/datenDownload.json"

r1 = session.get(url_wiski)
r2 = session.get(url_all)

data_wiski = r1.json()
data_all = r2.json()

https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/wiski/WISCHHAFEN_SPERRWERK!Salzgehalt.zip

# Filter variable
variable='Salzgehalt.zip'
Salzgehalt=[entry for entry in data_all if variable in entry['zipFileName']]

zipfiles=[entry['zipFileName'] for entry in Salzgehalt]

for zipfile in zipfiles:
    url='/'.join((base,zipfile))

try:
    r = requests.head(url, allow_redirects=True, timeout=5)
    if r.status_code == 200:
        print(" ✅ exists — downloading", os.path.basename(url))
        data = requests.get(url)
        with open(os.path.basename(url), "wb") as f:
            f.write(data.content)
    else:
        print(" ❌ not found:", r.status_code)
except Exception as e:
    print(" ⚠️ error:", e)

time.sleep(0.5)





# 2. Fetch mapping CSV for parameters
url_map = "https://www.kuestendaten.de/media/zdm/WISKI/Mapping_Parametergruppe_Parameter.csv"
r3 = session.get(url_map)
r3.raise_for_status()
decoded = r3.content.decode('utf-8').splitlines()
param_mapping = list(csv.reader(decoded, delimiter=';'))

# 3. Inspect the JSON structures (you may print a few entries)
print(data_wiski[:2])
print(data_all[:2])

# 4. Write logic to build URLs
results = []
for rec in data_all:
    station = rec.get("Stationenname", rec.get("stationName"))  # adjust keys based on JSON
    folder = rec.get("Folder", "")  # example
    linkId = rec.get("DownloadID", rec.get("linkId"))
    parameters = rec.get("Parameter", [])
    for param in parameters:
        url = f"{base}/{folder}/{linkId}!{param}.zip" if folder else f"{base}/{linkId}!{param}.zip"
        results.append({"station": station, "linkId": linkId, "parameter": param, "url": url})

# 5. Save results
with open("download_links.csv", "w", newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=["station","linkId","parameter","url"])
    writer.writeheader()
    writer.writerows(results)

# 6. Optional: download each & log status
for rec in results:
    try:
        r = session.head(rec['url'], allow_redirects=True, timeout=10)
        rec['status'] = r.status_code
    except Exception as e:
        rec['status'] = f"error: {e}"
    time.sleep(0.2)

# Optionally write status column too.
