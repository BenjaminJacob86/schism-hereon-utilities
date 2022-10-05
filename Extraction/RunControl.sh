#!/bin/bash


# Settings
rundir=/work/gg0028/g260114/RUNS/Europe/europe/baroclinic
schism=/work/gg0028/g260114/RUNS/RunControl/codes/pschism_DKRZ_mpiifor_intel18_intelmpi_O2_EVAP_SB
batchfile=/work/gg0028/g260114/RUNS/RunControl/subroutines/run_schism
hotcombine=/work/gg0028/g260114/svn/trunk/src/Utility/Combining_Scripts/combine_hotstart7neu.XX
speedup_py=/work/gg0028/g260114/RUNS/RunControl/subroutines/speedup.py

nproc=1080
ntracer=2
cdir=$rundir/outputs
maindir=$(pwd)

 # computation
 nodes=30	
 taskspernode=36
 partition=compute2	 


createhot=1

############### S T A R T ####################
runfin=0

# ,ake check for ihot

# set hpc settings #############################################

#SBATCH --partition=comppart   # Specify partition name
#SBATCH --nodes=nnodes         # Specify number of nodes as 
#SBATCH --ntasks-per-node=ntask  # Specify number of tasks on each node
#SBATCH --time=08:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --mail-user=benjamin.jacob@hzg.de  # Set your e-mail address
#SBATCH --account=gg0028       # Charge resources on this project account
#SBATCH --output=pschism.o%j       # File name for standard output
#SBATCH --error=pschism.e%j        # File name for standard error output

cd subroutines
cp run_schism.sh0 run_schism
sed -i "s/partition=.*#/partition=$partition #/" run_schism
sed -i "s/nodes=.*#/nodes=$nodes #/" run_schism
sed -i "s/ntasks-per-node=.*#/ntasks-per-node=$taskspernode #/" run_schism
cd ..
########################################################################


################## R U N    M O D E L ################################# 
# repeat   [make hotstart] [run model] while run not cpmpleted
count=0
#while ((runfin==0)); do { 
while ((count<12)); do {
((count++))
############## C H E C K  R U N C O M P L E T E ################### 

cd $rundir
FILE=mirror.out
if [ -f $FILE ]; then

# test run for completion
if grep  "Run completed successfully" mirror.out
then
    # code if found
        echo run completed
	runfin=1
        exit
else
    # code if not found
    echo run not completed
fi
else
    # code if not found
    echo run not completed
fi
cd $maindir



if (( createhot == 1 ))
 then {
######   C R E A T E    L A S T   H o T S T A R T #######
# get comandline arguments from calling script:
# 1) schism setup directory
# 2) source code: combining script path
# 3) it: iterationstep
# 4) nproc: nuber of processors
# 5) ntracer: number of tracers used (e.g. for salt and temp only)
if (( count == 1 )) 
then {
RES=$(sbatch subroutines/createLastHot.sh $rundir $hotcombine $nproc $ntracer )
myJobID=${RES##* }
}
else {
RES=$(sbatch --dependency=afterok:$runJobID subroutines/createLastHot.sh $rundir $hotcombine $nproc $ntracer )
myJobID=${RES##* }
}
fi
echo job id=$myJobID 
##########################################################
}
fi
createhot=1



######### run model after hotstart combined
 # restart Model with hotstart if creation succeded

# set hpc settings

if (( createhot == 0 ))
 then {

RES=$(sbatch $batchfile)
runJobID=${RES##* }

}
else {
RES=$(sbatch --dependency=afterok:$myJobID $batchfile)
runJobID=${RES##* }
}
fi



# get latest combined stack
fname=$(find -name 'schout_0000*nc' | sort -V | tail -1) #get file with highest hotsart iteration number 
fname=${fname:2}                               # remove "./" from fname
stack1=$(echo ${fname} | sed 's/.*_//'  | sed 's/\.*.nc//')
echo last stackoutput is $stack1 >>../infotext


#------------------------

 }; done
