

https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/wiski/ABBENFLETH_SPERRWERK!Salzgehalt.zip
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/wiski/ACHTHOEFENER_FLETH_SIEL!Salzgehalt.zip
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/belum!Salzgehalt.zip
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/wiski/BELUM_PEGEL!Salzgehalt.zip



import requests
from bs4 import BeautifulSoup
import os

page_url = "https://www.kuestendaten.de/DE/Services/Messreihen_Dateien_Download/Download_Zeitreihen_node.html"

resp = requests.get(page_url)
resp.raise_for_status()

soup = BeautifulSoup(resp.text, "html.parser")

links = []
for a in soup.find_all("a", href=True):
    href = a["href"]
    # check for ZIP files or the file-pattern you expect
    if href.lower().endswith(".zip"):
        # some links may be relative
        if href.startswith("/"):
            full = "https://www.kuestendaten.de" + href
        elif href.startswith("http"):
            full = href
        else:
            full = os.path.join(os.path.dirname(page_url), href)
        links.append(full)

print(f"Found {len(links)} ZIP links")

for link in links:
    filename = os.path.basename(link)
    print("Filename:", filename, " URL:", link)

# --- if you want to download ---
download = True
if download:
    for link in links:
        fname = os.path.basename(link)
        print("Downloading", fname)
        r = requests.get(link)
        r.raise_for_status()
        with open(fname, "wb") as f:
            f.write(r.content)
    print("All done.")
