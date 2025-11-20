import requests
from bs4 import BeautifulSoup
import os
import requests, os, time
import unicodedata


base = "https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download"
parameter = "Salzgehalt"
elbe_stations = [
    "Abbenfleth Sperrwerk",	
    "Achthöfener Fleth Siel",	
    "Belum",	
    "Belum Pegel",	
    "BHV Alter Leuchtturm",
    "Blankenese Anleger (WGMN)",
    "Blexen",	
    "Brake",	
    "Brunsbüttel (MPM)", 
    "Cuxhaven Alte Liebe",
    # Dastations
    "D1 - Hanskalbsand, Oberfläche, Elbe-km 643,0",	
    "D2 - Juelsand, Sohle, Elbe-Km 651,3",
    "D3 - Pagensand-Nord, Oberfläche, Elbe-Km 664,7",
    "D3 - Pagensand-Nord, Sohle, Elbe-Km 664,7",
    "D4 - Rhinplate-Nord, Oberfläche, Elbe-Km 676,5",
    "D4 - Rhinplate-Nord, Sohle, Elbe-Km 676,5",
# add more station names    
]

https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/bhv_brake!Salzgehalt.zip

def reformat_name(stat):
    # Remove/replace problematic characters
    s = stat.strip()
    s = s.replace('ü','ue').replace('ö','oe').replace('ä','ae').replace('ß','ss')
    # remove parentheses and contents
    s = s.replace('(','').replace(')','')
    # replace spaces and slashes
    s = s.replace(' ','_').replace('/','_').replace('\\','_')
    # uppercase everything
    s = s.upper()
    return s

stations = [reformat_name(stat) for stat in elbe_stations]


def clean(s):
    # lowercase and strip
    s = s.lower().strip()

    # remove number prefixes like "d1 -"
    if s.startswith("d") and "-" in s.split()[0]:
        s = s.split("-", 1)[1].strip()

    # replace umlauts
    repl = {
        "ä":"ae","ö":"oe","ü":"ue","ß":"ss"
    }
    for k,v in repl.items():
        s = s.replace(k, v)

    # remove accents just in case
    s = ''.join(c for c in unicodedata.normalize('NFKD', s) if not unicodedata.combining(c))

    # remove punctuation
    for char in ",()":
        s = s.replace(char, "")

    # space → underscore
    s = s.replace(" ", "_")

    return s


for rs in elbe_stations:
    key = rs.split()[0]

    if key.startswith("D"):  # D-profile
        name = clean(rs)
        name = "LZ_" + name
        url = f"{base}/{name}!{parameter}.zip"
    else:  # normal station
        name = clean(rs).upper()
        url = f"{base}/wiski/{name}!{parameter}.zip"

    print("Testing:", url)

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




https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/LZ_hanskalbsand_sohle!Salzgehalt.zip



https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/belum!Salzgehalt.zip


name="D1 - Hanskalbsand, Oberfläche, Elbe-km 643,0"
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/LZ_hanskalbsand_oberflaeche!Salzgehalt.zip

D1 - Hanskalbsand, Oberfläche, Elbe-km 643,0	
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/LZ_hanskalbsand_oberflaeche!Salzgehalt.zip
D2 - Juelsand, Oberfläche, Elbe-Km 651,3	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
https://www.kuestendaten.de/DE/dynamisch/appl/data/daten_prod/prodNeuProd/direct_download/LZ_juelsand_oberflaeche!Salzgehalt.zip



Dornbusch Brücke	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Dwarsgat	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Elbe 1 Sohle (BSH)	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Elbe 1 Wassersäule (BSH)	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Elsfleth	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Este Inneres Sperrwerk AP	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Farge	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Freiburg Sperrwerk	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Geversdorf Brücke	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Hahnöfer Sand West Siel	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
HB-HEMELINGEN (LF)	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Helgoland Ost Sohle (BSH)	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Helgoland Ost Wassersäule (BSH)	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Helgoland Süd Sohle (BSH)	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Helgoland Süd Wassersäule (BSH)	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Krückau-Sperrwerk BP	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Leuchtturm Alte Weser	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Lüheort	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Lühort	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
LZ 1b Krummendeich, Elbe-km 697,0	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
LZ 2a Neufeldreede	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
LZ 3c - Altenbrucher Bogen	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
LZ 4c Spitzsand West	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Neuenseer Schleusenfleth Siel	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Nordenham	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Osteriff (MPM) Sohle	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Otterndorf (MPM) Sohle	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Pinnau-Sperrwerk BP	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Rechtenfleth	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Robbensüdsteert	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Ruthenstrom-Sperrwerk	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Schöneworth Siel	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Schwinge-Sperrwerk	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Stör-Sperrwerk BP	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Strohhauser Plate Ost	Wasserstraßen- und Schifffahrtsamt Weser-Jade-Nordsee	+
Twielenfleth Siel	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	+
Wischhafen Sperrwerk	Wasserstraßen- und Schifffahrtsamt Elbe-Nordsee	-