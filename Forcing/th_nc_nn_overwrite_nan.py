# nan_fill.py
# extrapolate boundary data with nan (i.e. due to coarse forcing model) by overweriting nan with nearest neighbour vlaid values within th.nc

import glob
import xarray as xr
from matplotlib import pyplot as plt

thncs=glob.glob('*.th.nc')

for file in thncs:
	file
	break
