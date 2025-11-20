expname=ICM   #ClimateProj  #CNTRL #_no_wave #Veg_REF  #dtry10cm #oldcode #HydroOnly # Vegetation_Max  #dt90 #AaronImp #test_hot_control  #FRCww3BetaMax1.8
hotstart0=hotstart.nc_sed_merged #hotstart.nc
#bdspec=www.spec.nc

mkdir $expname
cd $expname
mkdir outputs hotstarts #combined
cp ../*.nml .
cp ../*.in .
cp ../*.sh .
cp ../*batch* .
cp -r scripts .

ln -s ../$hotstart0 hotstart.nc
ln -s ../*th.nc .
ln -s ../*nu.nc .
ln -s ../*.gr3 .
ln -s ../*grid* .
ln -s ../*.th* .
ln -s ../*.prop .
ln -s ../sflux .
ln -s ../*.ic .
ln -s ../ww3.spec.nc .
