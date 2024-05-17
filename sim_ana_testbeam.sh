#!/bin/bash
#source $SHIP_CVMFS_SETUP_FILE
#source $FAIRSHIP_DIR/config.sh
#source $SHIP_CVMFS_SETUP_FILE
#source /cvmfs/ship.cern.ch/SHiP-2021/latest/setUp.sh
#source /afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/start_ali.sh
#cd /afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
# cd /afs/cern.ch/user/u/ursovsnd/SND
source start_ali.sh
# cd /afs/cern.ch/user/u/ursovsnd/cluster/cluster
#Condor_folder=/afs/cern.ch/user/e/edursov/private/SIMULATIONS/my_FairShip/muon/run_processing/
#cd "$Condor_folder"
#alienv enter FairShip/latest
set -ux
echo "Starting script."
RUN=$1
ProcId=$2
LSB_JOBINDEX=$((ProcId+1))
NTOTAL=$4
NJOBS=$3
NAME=$5


#echo $MUONS
#echo $DIR
#echo $SUB
N=$(( NTOTAL/NJOBS + ( LSB_JOBINDEX == NJOBS ? NTOTAL % NJOBS : 0 ) ))
FIRST=$(((NTOTAL/NJOBS)*(LSB_JOBINDEX-1)))

#python $SNDSW_ROOT/shipLHC/run_simSND.py  --Ntuple  -n $N -f /eos/experiment/sndlhc/MonteCarlo/FLUKA/muons_up/version1/unit30_Nm.root  --eMin 1.0 --output "$EOS_DATA"/"$DIR"/"$LSB_JOBINDEX"/
# python new_mufi_exp.py --Nstart $FIRST -N $N
python -m pip install pandas
python new_mufi_exp.py -r $RUN -p /eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/ -g geofile_full.Ntuple-TGeant4_nom.root --nStart $FIRST --nEvents $N
python $SNDSW_ROOT/shipLHC/run_digi -r $RUN -p /eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/ -g geofile_full.Ntuple-TGeant4_nom.root --nStart $FIRST --nEvents $N
# python new_mufi_exp.py -r $RUN -p /eos/user/e/ekhaliko/Documents/SND_Data/test_100GeV_n10k_aug2023_pi+/ -f /eos/user/e/ekhaliko/Documents/SND_Data/test_100GeV_n10k_aug2023_pi+/merge_matei_cuts_no_thr_no_sat_wo-neutrals.root  -g geofile_full.PG_211-TGeant4.root --Nstart $FIRST --nEvents $N

# xrdcp $RUN root://eosuser.cern.ch//eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/"$RUN"_"$NAME"/$LSB_JOBINDEX/$RUN 	
xrdcp output.csv root://eosuser.cern.ch//eos/user/u/ursovsnd/SWAN_projects/tests/from_condor/"$RUN"_"$NAME"/$LSB_JOBINDEX/output.csv
# xrdcp *.root root://eosuser.cern.ch//eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/"$RUN"_"$NAME"/$LSB_JOBINDEX/output.root