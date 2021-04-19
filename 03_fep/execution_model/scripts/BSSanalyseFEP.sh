#!/usr/bin/env bash
echo $SLURM_JOB_NAME
date
#$BSSHOME/sire.app/bin/python scripts/BSSanalyseFEP.py $1 $2 $3

lig0=$1
lig1=$2
engine=$3

if [ $engine = "SOMD" ]
then
    echo "USING SOMD"
    #cd run/$lig0~$lig1/$stage/$lambda/
    echo "cd run/$lig0~$lig1/bound"
    #$BSSHOME/sire.app/bin/analyse-freenrg mbar -i */simfile.dat -p 95 --overlap --subsampling -o freenrg.dat 1> analyse_free_nrg.log 2> analyse_free_nrg.err
    echo "$BSSHOME/sire.app/bin/analyse-freenrg mbar -i */simfile.dat -p 95 --overlap --subsampling -o freenrg.dat 1> analyse_free_nrg.log 2> analyse_free_nrg.err"
    #cd ../free
    echo "cd ../free"
    #$BSSHOME/sire.app/bin/analyse-freenrg mbar -i */simfile.dat -p 95 --overlap --subsampling -o freenrg.dat 1> analyse_free_nrg.log 2> analyse_free_nrg.err
    echo "$BSSHOME/sire.app/bin/analyse-freenrg mbar -i */simfile.dat -p 95 --overlap --subsampling -o freenrg.dat 1> analyse_free_nrg.log 2> analyse_free_nrg.err"
    # Extract free energy from bound/freenrg.dat and free/freenrg.dat and compute relative binding free energy
    #cd ../
    echo "cd ../"
    # TODO with dummy MBAR outputs
    bindingnrg=-0.40 # kcal/mol
    bindingnrgerr=0.20 # kcal/mol
elif [ $engine = "GROMACS" ]
then
    echo "USING GROMACS"
else
    echo "The FEP engine $engine is not supported"
fi

# Deposit relative binding free energy in output folder (check if not exists already?)
echo "$lig0~$lig1,$bindingnrg,$bindingnrgerr" >> output/freeenergies.csv


sleep 5
exit 0
