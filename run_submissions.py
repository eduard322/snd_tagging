import os
import sys


#SUBFOLDER = os.environ["PRODUCTION_FOLDER"]
#DATA_DIR = os.path.join(os.environ["EOS_DIR"], SUBFOLDER)
N_JOBS= 100
DIR = sys.argv[1]
n_events = 100000
eos = os.environ['EOS_DATA']
#with open("input_for_muon_prod.txt", "r") as f:
'''
    for line in f:
        filepath, n_events, foldername = line.strip().split(", ")
        directory = os.path.join(DATA_DIR, foldername)
        if not os.path.exists(directory): 
            os.makedirs(directory)
'''

'''
for i in range(N_JOBS):
    directory = os.path.join(DIR, str(i+1))
    if not os.path.exists(directory):
	os.makedirs(directory)
'''

	
os.system("condor_submit directory={dir} N={n_jobs} n_events={ne} sim.sub".format(dir=DIR, n_jobs=N_JOBS, ne=n_events))
#os.system("hadd {dir}/total.root *_rec.root".format(dir=eos))
#os.system("python new_ana_up.py -f {dir}/total.root -g {dir}/geofile_full.conical.Pythia8-TGeant4.root".format(dir=eos))
