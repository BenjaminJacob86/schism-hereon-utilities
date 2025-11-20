import requests, os, time

session = requests.Session()
base = "https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download"

# Load WISKI metadata
url = f"{base}/wiski/WISKI_datenDownload.json"
data = session.get(url).json()


# 1. Fetch JSONs
url_wiski = base + "/wiski/WISKI_datenDownload.json"
url_all = base + "/datenDownload.json"

r1 = session.get(url_wiski)
r2 = session.get(url_all)

data_wiski = r1.json()
data_all = r2.json()


https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/belum!Salzgehalt.zip
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/wiski/BELUM_PEGEL!Salzgehalt.zip

# Filter variable 
variable='Salzgehalt.zip'
Salzgehalt1=[entry for entry in data_wiski if variable in entry['zipFileName']]
Salzgehalt2=[entry for entry in data_all if variable in entry['zipFileName']]


outdir = "downloads"
os.makedirs(outdir, exist_ok=True)

def try_download(fname):
    # try paths in order
    paths = [
        f"{base}/wiski/{fname}",
        f"{base}/{fname}"
    ]

    for url in paths:
        print(f"🔎 Checking {url}")
        head = session.head(url, allow_redirects=True, timeout=10)

        if head.status_code == 200:
            print(f"✅ Found: {url}")
            r = session.get(url, stream=True)

            local = os.path.join(outdir, fname.replace("!","_"))
            with open(local, "wb") as f:
                for chunk in r.iter_content(8192):
                    f.write(chunk)

            print(f"📦 Saved: {local}")
            return True
        else:
            print(f"❌ {url} -> {head.status_code}")

    print(f"⚠️ Failed all paths for {fname}")
    return False


for rec in Salzgehalt: #data:
    fname = rec["zipFileName"]
    try_download(fname)
    time.sleep(0.2)


for rec in Salzgehalt2: #data:
    fname = rec["zipFileName"]
    print(fname)
    try_download(fname)
    time.sleep(0.6)


test=[item for item in Salzgehalt2 if 'heme' in item['zipFileName']]
for rec in test: #data:
    fname = rec["zipFileName"]
    print(fname)
    try_download(fname)
    time.sleep(0.6)

https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/Helgoland_Ost_Wasser!Salzgehalt.zip



https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/bhv_rechtenfleth!Salzgehalt.zip
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/bhv_rechtenfleth!Salzgehalt.zip


https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/bhv_lt_hb_hemelingen!Salzgehalt.zip 
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/bhv_lt_hb_hemelingen!Salzgehalt.zip

for rec in Salzgehalt: #data:
    print(rec["zipFileName"])
    
    
    
    
   paths = [
        f"{base}/{fname}"
    ]

    for url in paths:
        print(f"🔎 Checking {url}")
        head = session.head(url, allow_redirects=True, timeout=10)
        
        
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/wiski/WISCHHAFEN_SPERRWERK!Salzgehalt.zip