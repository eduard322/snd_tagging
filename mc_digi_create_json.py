import json
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument(
#    'inputfile',
#    help='''Simulation results to use as input. '''
#    '''Supports retrieving files from EOS via the XRootD protocol.''')

parser.add_argument(
    '--geofile',
    dest="geo",
    type = str,
    default="geo.root")  

   
args = parser.parse_args()

# Extract all arguments as a dictionary
arguments_dict = vars(args)
json_file_path = f'arguments_{args.koef:.2f}.json'
with open(f"params/{json_file_path}", 'w') as json_file:
    json.dump(arguments_dict, json_file, indent=4)
