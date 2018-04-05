import pandas as pd
import os
import re
import argparse
import networkx as nx
import sys
import pickle

parser = argparse.ArgumentParser('Reads .mol2 files and returns Orca 4.0 input files.')
parser.add_argument('-i', '--input_mol2', help='input file name')
parser.add_argument('-o', '--output', default='output', help='The desired name of output files.')
parser.add_argument('-ip', '--ipath', default=os.getcwd(), help='Directory of the .mol2 file.')
parser.add_argument('-op', '--opath', default=os.getcwd(), help='Directory for the output files.')
parser.add_argument('-m', '--molecules', default=2, help='The number of different molecules in your system.')

args = parser.parse_args()

in_path = args.ipath
out_path = args.opath
name = args.input_mol2
in_file = os.path.join(in_path, name)

file_format = re.search('\.mol2', name)
if file_format is None:
    print('input not a mol2 file!!!')
    sys.exit(1)

out_name = args.output
out_ext = {'orca': '.inp'}
os.makedirs(os.path.dirname(args.opath), exist_ok=True)

mol_types = {}
for i in range(0, args.molecules, 1):
    key = 'mol'+str(i)
    print('How many atoms are there in molecule '+str(i))
    while True:
        mol_atoms = input(':')
        try:
            mol_types[key] = int(mol_atoms)
            break
        except ValueError:
            print('Please enter an integer.')
        continue

with open(in_file, mode='r') as f:
    count = 0
    at_start,meta_start, bd_start = 0, 0, 0
    for line in f:
        if re.search('(.*)MOLECULE(.*)', line):
            meta_start = count
        elif re.search('(.*)ATOM(.*)', line):
            at_start = count
        elif re.search('(.*)BOND(.*)', line):
            bd_start = count
        count += 1

meta = pd.read_csv(in_file,
                    header=None,
                    sep=' ',
                    skiprows=meta_start+2,
                    nrows=1
                    )
with open('meta', mode='wb') as f:
    pickle.dump(meta, f)
    
atoms = pd.read_csv(in_file, 
                    header=None, 
                    sep=' ', 
                    skiprows=at_start+1, 
                    nrows=meta.loc[0, 0]
                    )
with open('atoms', mode='wb') as f:
    pickle.dump(atoms, f)

bonds = pd.read_csv(in_file, 
                    header=None,
                    sep=' ',
                    skiprows=bd_start+1,
                    nrows=meta.loc[0, 1]
                    )
with open('bonds', mode='wb') as f:
    pickle.dump(bonds, f)
