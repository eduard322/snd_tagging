#!/bin/bash
#source $SHIP_CVMFS_SETUP_FILE
#source $FAIRSHIP_DIR/config.sh
#source $SHIP_CVMFS_SETUP_FILE
source /cvmfs/sndlhc.cern.ch/latest/setUp.sh
Condor_folder=/afs/cern.ch/user/u/ursovsnd/SND
source $Condor_folder/start_ali.sh

#alienv enter FairShip/latest
set -ux
echo "Starting script."
DIR=$1
ProcId=$2
LSB_JOBINDEX=$((ProcId+1))
NTOTAL=$4
NJOBS=$3
N=$(( NTOTAL/NJOBS + ( LSB_JOBINDEX == NJOBS ? NTOTAL % NJOBS : 0 ) ))
FIRST=$(((NTOTAL/NJOBS)*(LSB_JOBINDEX-1)))
if eos stat "$EOS_DATA"/"$DIR"/"$LSB_JOBINDEX"/ship.conical.MuonBack-TGeant4.root; then
	echo "Target exists, nothing to do."
	exit 0
else
	
	#python "$FAIRSHIP"/macro/run_simScript.py --muShieldDesign $MUSHIELD --MuonBack --nEvents $N --firstEvent $FIRST -f $MUONS --FastMuon -g remove_0_X.root --stepMuonShield --output "$EOS_DATA"/"$DIR"/"$SUB"/"$LSB_JOBINDEX"/
	#xrdcp ship.conical.MuonBack-TGeant4.root root://eospublic.cern.ch/"$EOS_DATA"/"$DIR"/"$LSB_JOBINDEX"/ship.conical.MuonBack-TGeant4.root
	#python FLUX_4.py "$EOS_DATA"/"$DIR"/"$SUB"/"$LSB_JOBINDEX"/ship.conical.MuonBack-TGeant4.root -n "$EOS_DATA"/"$DIR"/"$SUB"/"$LSB_JOBINDEX"  
	#python flux_map_custom.py "$EOS_DATA"/"$DIR"/"$SUB"/"$LSB_JOBINDEX"/ship.conical.MuonBack-TGeant4.root -n "$EOS_DATA"/"$DIR"/"$SUB"/"$LSB_JOBINDEX" -o "$EOS_DATA"/"$DIR"/"$SUB"/"$LSB_JOBINDEX"/flux_map.root   
	python /afs/cern.ch/user/u/ursovsnd/cluster/cluster/flux_map_custom.py --Eloss "$ELOSS" --Pcut "$PCUT"
	xrdcp flux_map.root root://eosuser.cern.ch//eos/user/e/ekhaliko/Documents/SND_Data/test_E100-2500_n30k_center_pi-/"$DIR"/"$LSB_JOBINDEX"/flux_map_"$PCUT"_"$ELOSS".root 	

fi
