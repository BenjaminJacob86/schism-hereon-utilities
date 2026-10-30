"""
EasyGSH → SCHISM sediment mapping workflow

Steps:
1. Read EasyGSH raster sediment distributions (phi space)
2. Compute class fractions (Pi)
3. Aggregate into SCHISM sediment classes
4. Rasterize structured grid (100 m)
5. Bilinear interpolation to SCHISM mesh
6. Depth-based extrapolation for missing areas
7. Normalize + export bed fractions


Used dataset: 
https://mdi-de.baw.de/easygsh/Easy_DownloadS.html#home
2016
Titel	Format	Vorschau	Metadaten	Download
Hauptkomponenten (Kurzform), Deutsche Bucht	Shape	icon auge schwarz	Datenrepository	ZIP
Hauptkomponenten (Langform), Deutsche Bucht	Shape	icon auge schwarz	Datenrepository	ZIP
Nebenkomponenten (Kurzform), Deutsche Bucht	Shape	icon auge schwarz	Datenrepository	ZIP
Nebenkomponenten (Langform), Deutsche Bucht	Shape	icon auge schwarz	Datenrepository	ZIP
Titel	Format	Vorschau	Metadaten	Download
Median-Korndurchmesser d50, Deutsche Bucht, 100m Raster	GeoTIFF	icon auge schwarz	Datenrepository	ZIP
Median-Korndurchmesser d50 Isolinien, Deutsche Bucht	Shape	icon auge schwarz	Datenrepository	ZIP
Sortierung, Deutsche Bucht, 100m Raster	GeoTIFF	icon auge schwarz	Datenrepository	ZIP
Schiefe, Deutsche Bucht, 100m Raster	GeoTIFF	icon auge schwarz	Datenrepository	ZIP
Porosität, Deutsche Bucht, 100m Raster	GeoTIFF	icon auge schwarz	Datenrepository	ZIP
Titel	Format	Vorschau	Metadaten	Download
Deutsche Bucht	CSV	-	Datenrepository	ZIP
2016 - komplett	-	-	Datenrepository	ZIP

Should be in the CSV files used file:
'100mRasterDB_2021.csv

"""

# =========================================================
# 0. IMPORTS & SETUP
# =========================================================
import numpy as np
import pandas as pd
import os
from scipy.interpolate import griddata

# Add custom SCHISM library  functions
# https://github.com/BenjaminJacob86/schism-hereon-utilities
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code//schism/scripts/schism-hereon-utilities/')
from schism import *

####

# =========================================================
# 0. SETTINGS
# =========================================================
schism_dir='/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_sed_routine_seagrass/'

# =========================================================
# 1. READ EASYGSH DATA
# =========================================================
Sed = pd.read_csv('/gpfs/work/ksddata/EasyGsh/100mRasterDB_2021.csv', sep=';', skipinitialspace=True)

x = pd.to_numeric(Sed.iloc[:, 3], errors='coerce')
y = pd.to_numeric(Sed.iloc[:, 4], errors='coerce')
z = pd.to_numeric(Sed.iloc[:, 5], errors='coerce')
d50_given = pd.to_numeric(Sed.iloc[:, 7], errors='coerce')

# Cumulative grain size distribution (CDF)
F = Sed.iloc[:, 8:].apply(pd.to_numeric, errors='coerce').to_numpy()

# φ axis
phi_values = Sed.columns[8:].astype(float).to_numpy()
d_values = 2 ** (-phi_values)

# =========================================================
# 2. DERIVE CLASS FRACTIONS (PDF FROM cumulative CDF)
# =========================================================
Pi = -(np.diff(F, axis=1))   # shape: (n_points, n_classes)

# =========================================================
# 3. VALIDATION: d50 CHECK
# =========================================================
idx = np.argmax(F < 0.5, axis=1)

phi50 = (
    phi_values[idx - 1]
    + (0.5 - F[np.arange(len(F)), idx - 1])
    * (phi_values[idx] - phi_values[idx - 1])
    / (F[np.arange(len(F)), idx] - F[np.arange(len(F)), idx - 1])
)

d50_calc = 2 ** (-phi50) # should approximately the delived d50 (d50_given)

# =========================================================
# 4. READ SCHISM SEDIMENT CLASSES
# =========================================================
os.chdir('/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_sed_routine_seagrass/NmL_Johannes')

with open('sediment.nml') as f:
    line = ''
    while 'Sd50' not in line:
        line = f.readline()

d50_schism = np.asarray(
    ''.join(line.replace('\n','').replace('d','e').split('=')[1:]).split(','),
    np.float32
)

phi_schism = -np.log2(d50_schism)

# =========================================================
# 5. BUILD PHI BIN EDGES (FINE GRID)
# =========================================================
dphi = np.diff(phi_values)[0]

phi_edges = np.concatenate([
    [phi_values[0] - dphi/2],
    phi_values[:-1] + dphi/2
])

# =========================================================
# 6. SCHISM CLASS EDGES
# =========================================================
phi_s = np.sort(phi_schism)

phi_s_edges = np.concatenate([
    [-np.inf],
    (phi_s[:-1] + phi_s[1:]) / 2,
    [np.inf]
])

# =========================================================
# 7. AGGREGATE TO SCHISM CLASSES
# =========================================================
n_points = Pi.shape[0]
n_coarse = len(phi_s_edges) - 1

Pi_coarse = np.zeros((n_points, n_coarse))

for i in range(Pi.shape[1]):
    fine_left = phi_edges[i]
    fine_right = phi_edges[i + 1]
    fine_width = fine_right - fine_left

    for j in range(n_coarse):
        coarse_left = phi_s_edges[j]
        coarse_right = phi_s_edges[j + 1]

        overlap = max(
            0,
            min(fine_right, coarse_right) - max(fine_left, coarse_left)
        )

        if overlap > 0:
            Pi_coarse[:, j] += Pi[:, i] * (overlap / fine_width)

# sanity check
assert np.allclose(Pi_coarse.sum(axis=1), 1.0, equal_nan=True)

# =========================================================
# 8. STRUCTURED GRID (RASTERIZATION)
# =========================================================
x_unique = np.unique(x)
y_unique = np.unique(y)

Nx, Ny = len(x_unique), len(y_unique)

x_idx = np.searchsorted(x_unique, x)
y_idx = np.searchsorted(y_unique, y)

# build structured grid containers
grid = np.full((Ny, Nx), np.nan)
Pi_grid = np.full((Ny, Nx, n_coarse), np.nan)

# fill grid
grid[y_idx, x_idx] = d50_given
Pi_grid[y_idx, x_idx, :] = Pi_coarse

# =========================================================
# 9. INTERPOLATE TO SCHISM MESH (BILINEAR)
# =========================================================
os.chdir(schism_dir)
s = schism_setup()

X_grid, Y_grid = np.meshgrid(x_unique, y_unique)

mask = ~np.isnan(Pi_grid[:, :, 0])

xv = X_grid[mask]
yv = Y_grid[mask]

Pi_interp = np.full((len(s.x), n_coarse), np.nan)

for k in range(n_coarse):
    values = Pi_grid[:, :, k][mask]

    Pi_interp[:, k] = griddata(
        (xv, yv),
        values,
        (s.x, s.y),
        method='linear'
    )

# =========================================================
# 10. DEPTH-BASED EXTRAPOLATION (ESTUARIES / GAPS)
# =========================================================
depth = np.asarray(s.depths)

missing = np.isnan(Pi_interp).any(axis=1)

bins = [-np.inf,0,5,10,15,20,25,30,35,40,45,50,np.inf]

Pi_filled = Pi_interp.copy()

for k in range(n_coarse):
    for i in range(len(bins) - 1):

        dmin, dmax = bins[i], bins[i + 1]

        sel_valid = (~missing) & (depth > dmin) & (depth <= dmax)
        sel_missing = (missing) & (depth > dmin) & (depth <= dmax)

        if np.any(sel_valid):
            mean_val = np.nanmean(Pi_interp[sel_valid, k])
            Pi_filled[sel_missing, k] = mean_val

# =========================================================
# 11. FINAL FALLBACK (GLOBAL MEAN)
# =========================================================
still_nan = np.isnan(Pi_filled).any(axis=1)

if np.any(still_nan):
    global_mean = np.nanmean(Pi_interp, axis=0)
    Pi_filled[still_nan] = global_mean

# =========================================================
# 12. NORMALIZATION (CRITICAL)
# =========================================================
Pi_filled /= Pi_filled.sum(axis=1, keepdims=True)

# SCHISM ordering: inverse phi order
Pi_schism = Pi_filled[:, ::-1]

# =========================================================
# 13. EXPORT TO SCHISM
# =========================================================
for k in range(n_coarse):
    s.dump_gr3_spat_var(
        filename=f'bed_frac_{k+1}.ic',
        nodevalues=Pi_schism[:, k]
    )