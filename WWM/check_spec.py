import xarray as xr

ds1=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/WAVE/ww3.2017_cat_01_03.nc')
ds2=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/WAVE/FRCww3BetaMax1.8/ww3.2017cat0103.spec.nc')

ww3a=WW3_mesh(diri+ww3mesh)
ww3a.load_points(diri+ww3PointList)
ww3a.open_point_output('/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_from_mistral/BetaMax2.4/ww3.2017_cat_01_03.nc')

ww3b=WW3_mesh(diri+ww3mesh)
ww3b.load_points(diri+ww3PointList)
ww3b.open_point_output('/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_from_mistral/FRCWW31.8/ww3.2017cat0103.spec.nc')

#difference
nr=-1
plt.figure()
plt.plot(ww3a.Hs[:,nr])
plt.plot(ww3b.Hs[:,nr])
# bdforcing differ