import json
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument(
#    'inputfile',
#    help='''Simulation results to use as input. '''
#    '''Supports retrieving files from EOS via the XRootD protocol.''')

parser.add_argument(
    '--Energy',
    dest="Energy",
    type = int,
    default=100)  

parser.add_argument(
    '--Merged',
    dest="merged",
    type = bool,
    default=True)

parser.add_argument(
    '--Koef',
    dest="koef",
    type = float,
    default=True)     
args = parser.parse_args()

# Extract all arguments as a dictionary
arguments_dict = vars(args)
json_file_path = f'arguments_{args.koef:.2f}_{args.Energy}.json'
with open(f"params/{json_file_path}", 'w') as json_file:
    json.dump(arguments_dict, json_file, indent=4)
