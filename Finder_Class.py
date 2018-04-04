import pandas as pd
import os
import re
import argparse
import networkx as nx
import sys


class Parser:

    def __init__(self):

        self.args = arg_parser()

    def arg_parser(self):

        parser = argparse.ArgumentParser('Reads .mol2 files and returns Orca 4.0 input files.')
        parser.add_argument('-i', '--input_mol2', help='input file name')
        parser.add_argument('-o', '--output', default='output', help='The desired name of output files.')
        parser.add_argument('-ip', '--ipath', default=os.getcwd(), help='Directory of the .mol2 file.')
        parser.add_argument('-op', '--opath', default=os.getcwd(), help='Directory for the output files.')
        parser.add_argument('-m', '--molecules', default=2, help='The number of different molecules in your system.')

        args = parser.parse_args()

        return args


class Finder:

    def __init__(self, args):
        in_path = args.ipath
        out_path = args.opath
        # name = args.input_mol2
        name = 'Supercell_XYZ_NVT_20Ps.mol2'
        file_format = re.search('\.mol2', name)
        if file_format is None:
            print('input not a mol2 file!!!')
            sys.exit(1)

        self.in_file = os.path.join(in_path, name)



        out_name = args.output
        out_ext = {'orca': '.inp'}
        os.makedirs(os.path.dirname(args.opath), exist_ok=True)

        self.mol_types = {}
        for i in range(0, args.molecules, 1):
            key = 'mol' + str(i)
            print('How many atoms are there in molecule ' + str(i))
            while True:
                mol_atoms = input(':')
                try:
                    self.mol_types[key] = int(mol_atoms)
                    break
                except ValueError:
                    print('Please enter an integer.')
                    continue
