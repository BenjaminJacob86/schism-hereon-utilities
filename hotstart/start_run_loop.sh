#!/bin/bash
batchfile=run_levante_chain
# run schism model ntimes
# with creating hotstart inbetween

set +k
cdir=$(pwd)
cd ${cdir}
batchfile=run_levante_chain
# runnr == 0 means hotstart at time 0
count=$(<runnr)
for (( i=0;i<2;i++ ))
do
echo 
#if (( count==0 )) # first hotsart
if (( i==0 )) # first hotsart
then
	jid=$(sbatch $batchfile)
else	
	jid=$(sbatch --dependency=afterok:$jidHot $batchfile)
fi	
	jid=${jid:20}
	jidHot=$(sbatch --dependency=afterok:$jid $cdir/createLastHot.sh)
	jidHot=${jidHot:20}

done