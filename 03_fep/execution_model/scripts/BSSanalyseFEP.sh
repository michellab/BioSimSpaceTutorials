#!/usr/bin/env bash
echo $SLURM_JOB_NAME
date

lig0=$1
lig1=$2
engine=$3

cd ./outputs

if [[ $engine == *"SOMD"* ]]
then
    echo "USING SOMD"
    cd SOMD/$lig0~$lig1/bound/
    
    $BSSHOME/biosimspace.app/bin/analyse-freenrg mbar -i */simfile.dat -p 95 --overlap --subsampling -o freenrg.dat 1> analyse_free_nrg.log 2> analyse_free_nrg.err
   
    # Extract free energy from bound/freenrg.dat and compute relative binding free energy
    bindingnrg_bound=$(awk -F " " '/#MBAR free energy/{getline; print $1}' freenrg.dat | sed 's/,//g')
    bindingnrgerr_bound=$(awk -F " " '/#MBAR free energy/{getline; print $2}' freenrg.dat | sed 's/,//g')

    cd ../free
    
    $BSSHOME/biosimspace.app/bin/analyse-freenrg mbar -i */simfile.dat -p 95 --overlap --subsampling -o freenrg.dat 1> analyse_free_nrg.log 2> analyse_free_nrg.err
    
    # Extract free energy from free/freenrg.dat and compute relative binding free energy
    bindingnrg_free=$(awk -F " " '/#MBAR free energy/{getline; print $1}' freenrg.dat | sed 's/,//g')
    bindingnrgerr_free=$(awk -F " " '/#MBAR free energy/{getline; print $2}' freenrg.dat | sed 's/,//g')

    cd ../../
    

elif [ $engine = "GROMACS" ]
then
    echo "USING GROMACS"
    # code here to do the same as above, but in GROMACS.
else
    echo "The FEP engine $engine is not supported"
fi

# Deposit relative binding free energy in output/engine/. 
echo "$lig0~$lig1,$bindingnrg_bound,$bindingnrgerr_bound,$bindingnrg_free,$bindingnrgerr_free" >> freeenergies.csv


sleep 5
exit 0
