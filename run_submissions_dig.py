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
    '--param',
    dest="param",
    help='''run''',
    type=str)
parser.add_argument(
'-j',
'--jobs',
dest="njobs",
help='''number of jobs''',
type=int,
default=1)
parser.add_argument(
'-n',
'--nEvents',
dest="nEvents",
help='''number of events''',
type=int,
default=1)

parser.add_argument(
    '--outputGlobalFolder',
    default='/eos/user/u/ursovsnd/private/',
    dest='outputGlobalFolder',
    help='''File to write the flux maps to. '''
    '''Will be recreated if it already exists.''')
parser.add_argument(
    '--outputFolder',
    default='default',
    dest='outputFolder',
    help='''File to write the flux maps to. '''
    '''Will be recreated if it already exists.''')

    
args = parser.parse_args()

# input("you sure?")

if os.path.exists(f"./logs/{args.outputFolder}"): 
    shutil.rmtree(f"./logs/{args.outputFolder}")

os.makedirs(f"./logs/{args.outputFolder}")

outfold = os.path.join(args.outputGlobalFolder, args.outputFolder)

for i in range(args.njobs):
      os.makedirs(os.path.join(f"./logs/{args.outputFolder}", str(i)))
	
# if os.path.exists(f"/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/{args.run}_{args.name}/"):
#     shutil.rmtree(f"/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/{args.run}_{args.name}/")
# os.makedirs(f"/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/{args.run}_{args.name}/")
os.system("condor_submit arg_file={arg_file} N={n_jobs} nEvents={nEvents} out_fold={outfold} outputFolder={outputFolder} sub_files/dig.sub".format(arg_file=args.param, 
                                                                                                    n_jobs=args.njobs,
                                                                                                    nEvents=args.nEvents,
                                                                                                    outfold=outfold,
                                                                                                    outputFolder=args.outputFolder))
#os.system("hadd {dir}/total.root *_rec.root".format(dir=eos))
#os.system("python new_ana_up.py -f {dir}/total.root -g {dir}/geofile_full.conical.Pythia8-TGeant4.root".format(dir=eos))
