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
# name = 'Supercell_XYZ_NVT_20Ps.mol2'
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

with open(in_file, mode='r', encoding='utf-8') as f:
    met = 0
    at = 0
    bd = 0
    for line in f:
        met_test = re.search('MOLECULE', line)
        at_test = re.search('ATOM', line)
        bd_test = re.search('BOND', line)
        if met_test is not None:
            meta_start = met
        elif at_test is not None:
            at_start = at
        elif bd_test is not None:
            bd_start = bd
        else:
            continue
        met += 1
        at += 1
        bd += 1
        
    meta = pd.read_csv(in_file,
                      header=None,
                      delim_whitespace=True,
                      skiprows=meta_start+1,
                      nrows=1
                      )
    with open('meta', mode='wb') as f:
    	pickle.dump(meta, f)
    
    atoms = pd.read_csv(in_file, 
                        header=None, 
                        delim_whitespace=True, 
                        skiprows=at_start, 
                        nrows=meta.loc[0, 0]
                       )
    with open('atoms', mode='wb') as f:
    	pickle.dump(atoms, f)

    bonds = pd.read_csv(in_file, 
                        header=None,
                        delim_whitespace=True,
                        skiprows=bd_start,
                        nrows=meta.loc[0, 1]
                       )
    with open('bonds', mode='wb') as f:
    	pickle.dump(bonds, f)
