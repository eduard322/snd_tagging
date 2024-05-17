import os
import sys
import argparse

N_JOBS= 1

parser = argparse.ArgumentParser(description='Script to create flux maps.')
#parser.add_argument(
#    'inputfile',
#    help='''Simulation results to use as input. '''
#    '''Supports retrieving files from EOS via the XRootD protocol.''')
parser.add_argument(
    '-d',
    '--Dir',
    default='/eos/user/e/ekhaliko/Documents/SND_Data/test_E100-2500_n30k_center_pi-',
    help='''Directory''')
parser.add_argument(
'-P',
'--Pcut',
default= 0.0,
help='''set momentum cut''')
parser.add_argument(
'-E',
'--Eloss',
default= 0.0,
help= '''set Eloss cut''')

    
args = parser.parse_args()

	
os.system("condor_submit directory={dir} N={n_jobs} Eloss={Eloss} Pcut={Pcut} sim_ana.sub".format(dir=args.Dir, n_jobs=N_JOBS, Eloss=args.Eloss, Pcut = args.Pcut))
#os.system("hadd {dir}/total.root *_rec.root".format(dir=eos))
#os.system("python new_ana_up.py -f {dir}/total.root -g {dir}/geofile_full.conical.Pythia8-TGeant4.root".format(dir=eos))
