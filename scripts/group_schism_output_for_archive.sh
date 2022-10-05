#
max=371
window=45


i0=1
i1=$window

iwnd=1
while [ $iwnd -lt $((max/window)) ]
do echo $iwnd
outdir=jacobb_BS_2012to2017comb_part$iwnd
mkdir $outdir
for ((i=i0;i<=i1;i=i+1))
do 
	mv schout_$i.nc $outdir/
done 

iwnd=$((iwnd+1))
i0=$((i0+window))
i1=$((i1+window))
echo $i0
echo $i1
done


do echo $iwnd
i1=$max
echo $i0
echo $i1

outdir=jacobb_BS_2012to2017comb_part$iwnd
mkdir $outdir
for ((i=i0;i<=i1;i=i+1))
do 
	mv schout_$i.nc $outdir/
done 
