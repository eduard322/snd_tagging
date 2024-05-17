import os
import sys
import argparse
import shutil

parser = argparse.ArgumentParser(description='Script to create flux maps.')
#parser.add_argument(
#    'inputfile',
#    help='''Simulation results to use as input. '''
#    '''Supports retrieving files from EOS via the XRootD protocol.''')
parser.add_argument(
    '-r',
    '--run',
    dest="run",
    help='''run''',
    type=int)
parser.add_argument(
'-j',
'--jobs',
dest="njobs",
help='''number of jobs''',
type=int,
default=1)
parser.add_argument(
'-N',
'--nEvents',
dest="nEvents",
default=1,
help= '''number of events''',
type=int)
parser.add_argument(
'--name',
dest="name",
default=1,
help= '''name''')

    
args = parser.parse_args()

input("you sure?")

if os.path.exists("./logs/"): 
    shutil.rmtree("./logs/")

os.makedirs("./logs/")

for i in range(args.njobs):
      os.makedirs(os.path.join("logs", str(i)))
	
# if os.path.exists(f"/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/{args.run}_{args.name}/"):
#     shutil.rmtree(f"/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/{args.run}_{args.name}/")
# os.makedirs(f"/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/{args.run}_{args.name}/")
os.system("condor_submit run={dir} N={n_jobs} n_events={nEvents} name={name} sim_ana_testbeam.sub".format(dir=args.run, n_jobs=args.njobs, nEvents=args.nEvents, name=args.name))
#os.system("hadd {dir}/total.root *_rec.root".format(dir=eos))
#os.system("python new_ana_up.py -f {dir}/total.root -g {dir}/geofile_full.conical.Pythia8-TGeant4.root".format(dir=eos))
