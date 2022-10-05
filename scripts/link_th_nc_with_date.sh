#!/bin/bash
# link forcing files based on date name

date=20211101

rm ?????D.th.nc
rm uv3D.th.nc 


for file in ??????$date*.th.nc
do
ln -sf $file ${file:0:6}.th.nc
done


for file in uv3D$date*.th.nc
do
ln -sf $file ${file:0:4}.th.nc
done
