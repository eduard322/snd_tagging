#!/bin/bash
#source $SHIP_CVMFS_SETUP_FILE
#source $FAIRSHIP_DIR/config.sh
#source $SHIP_CVMFS_SETUP_FILE
#source /cvmfs/ship.cern.ch/SHiP-2021/latest/setUp.sh
#source /afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/start_ali.sh
#cd /afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/
source /cvmfs/sndlhc.cern.ch/latest/setUp.sh
cd /afs/cern.ch/user/u/ursovsnd/SND
source /afs/cern.ch/user/u/ursovsnd/SND/start_ali.sh
cd /afs/cern.ch/user/u/ursovsnd/cluster/cluster
#Condor_folder=/afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/muon/run_processing/
#cd "$Condor_folder"
#alienv enter FairShip/latest
set -ux
echo "Starting script."
DIR=$1
#SUB=$6
ProcId=$2
LSB_JOBINDEX=$((ProcId+1))
#MUONS=$4
NTOTAL=$4
NJOBS=$3
TANK=6
ISHIP=3
MUSHIELD=8

#echo $MUONS
#echo $DIR
#echo $SUB
N=$(( NTOTAL/NJOBS + ( LSB_JOBINDEX == NJOBS ? NTOTAL % NJOBS : 0 ) ))
FIRST=$(((NTOTAL/NJOBS)*(LSB_JOBINDEX-1)))

#python $SNDSW_ROOT/shipLHC/run_simSND.py  --Ntuple  -n $N -f /eos/experiment/sndlhc/MonteCarlo/FLUKA/muons_up/version1/unit30_Nm.root  --eMin 1.0 --output "$EOS_DATA"/"$DIR"/"$LSB_JOBINDEX"/
python sim_ana.py "$EOS_DATA"/"$DIR"/"$LSB_JOBINDEX"/sndLHC.Ntuple-TGeant4.root

