mkdir time_mismatch
for ((i=4;i<16;i++))
do
	a=$(printf '%02d' $i)
	mv wwm_out$a time_mismatch/ 
	mv outputs$a time_mismatch/
	mv hotstart.nc_$a time_mismatch/ 
	mv wwminput.nml_$a.nc time_mismatch/
	mv wwm_hotfile_out_$a.nc time_mismatch/
done
