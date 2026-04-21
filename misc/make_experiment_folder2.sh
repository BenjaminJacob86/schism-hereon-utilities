expname=VOR_TEST #Test_Seyed # RUN24f  #RUN_DWD2008  #CNTRL #_no_wave #Veg_REF  #dtry10cm #oldcode #HydroOnly # Vegetation_Max  #dt90 #AaronImp #test_hot_control  #FRCww3BetaMax1.8
hotstart0=hotstart.nc_sed_merged #hotstart.nc
#bdspec=www.spec.nc

mkdir $expname
cd $expname
#refdir=/work/gg0028/g260114/SETUPS/GhanaV3_2D/
refdir=/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/

mkdir outputs hotstarts #combined
cp $refdir/*.nml .
cp $refdir/*.in .
cp $refdir/*.sh .
cp $refdir/*batch* .
cp $refdir/*run* .
cp $refdir/*template* .


#ln -s ../$hotstart0 hotstart.nc
ln -s $refdir/*th.nc .
#ln -s $refdir/*nu.nc .
ln -s $refdir/*.gr3 .
ln -s $refdir/*grid* .
ln -s $refdir/*.th* .
ln -s $refdir/*.prop .
ln -s $refdir/sflux .
ln -s $refdir/*hot* .
#ln -s ../*.ic .
#ln -s ../ww3.spec.nc .

rm runnr
echo 00 > runnr

