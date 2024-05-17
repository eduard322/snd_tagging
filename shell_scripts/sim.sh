#!/bin/bash

#########################################################################
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
source start_ali.sh
set -ux
echo "Starting script."
ProcId=$1
NJOBS=$2
PARAM=$3
NTOTAL=$4
OUTFOLDER=$5
LSB_JOBINDEX=$((ProcId+1))

N=$(( NTOTAL/NJOBS + ( LSB_JOBINDEX == NJOBS ? NTOTAL % NJOBS : 0 ) ))
FIRST=$(((NTOTAL/NJOBS)*(LSB_JOBINDEX-1)))
#########################################################################

python -m pip install pandas
python *.py --nStart $FIRST --nEvents $N --param $PARAM
#python run_digiSND.py --nStart $FIRST --nEvents $N --param $PARAM
# if if test -f $OUTFOLDER/$LSB_JOBINDEX/output.csv; then
# 	echo "Target exists, REWRITING THE FILE."
#     rm root://eosuser.cern.ch/$OUTFOLDER/$LSB_JOBINDEX/output.csv   
#     xrdcp output.csv root://eosuser.cern.ch/$OUTFOLDER/$LSB_JOBINDEX/output.csv
# else
xrdcp output.csv root://eosuser.cern.ch/$OUTFOLDER/$LSB_JOBINDEX/output.csv
# fi