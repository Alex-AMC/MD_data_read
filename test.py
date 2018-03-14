import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def periodic(coords, vecs):
    # Centre
    initial = coords
    # Positions middle
    left = initial.copy()
    left['X'] = left['X'] - vecs.iloc[0, 0]
    left['Y'] = left['Y'] - vecs.iloc[0, 1]
    left['Z'] = left['Z'] - vecs.iloc[0, 2]

    right = initial.copy()
    right['X'] = right['X'] + vecs.iloc[0, 0]
    right['Y'] = right['Y'] + vecs.iloc[0, 1]
    right['Z'] = right['Z'] + vecs.iloc[0, 2]

    front = initial.copy()
    front['X'] = front['X'] - vecs.iloc[2, 0]
    front['Y'] = front['Y'] - vecs.iloc[2, 1]
    front['Z'] = front['Z'] - vecs.iloc[2, 2]

    lfront = front.copy()
    lfront['X'] = lfront['X'] - vecs.iloc[0, 0]
    lfront['Y'] = lfront['Y'] - vecs.iloc[0, 1]
    lfront['Z'] = lfront['Z'] - vecs.iloc[0, 2]

    rfront = front.copy()
    rfront['X'] = rfront['X'] + vecs.iloc[0, 0]
    rfront['Y'] = rfront['Y'] + vecs.iloc[0, 1]
    rfront['Z'] = rfront['Z'] + vecs.iloc[0, 2]

    back = initial.copy()
    back['X'] = back['X'] + vecs.iloc[2, 0]
    back['Y'] = back['Y'] + vecs.iloc[2, 1]
    back['Z'] = back['Z'] + vecs.iloc[2, 2]

    lback = back.copy()
    lback['X'] = lback['X'] - vecs.iloc[0, 0]
    lback['Y'] = lback['Y'] - vecs.iloc[0, 1]
    lback['Z'] = lback['Z'] - vecs.iloc[0, 2]

    rback = back.copy()
    rback['X'] = rback['X'] + vecs.iloc[0, 0]
    rback['Y'] = rback['Y'] + vecs.iloc[0, 1]
    rback['Z'] = rback['Z'] + vecs.iloc[0, 2]

    # Top row
    up = initial.copy()
    up['X'] = up['X'] + vecs.iloc[1, 0]
    up['Y'] = up['Y'] + vecs.iloc[1, 1]
    up['Z'] = up['Z'] + vecs.iloc[1, 2]

    upleft = up.copy()
    upleft['X'] = upleft['X'] - vecs.iloc[0, 0]
    upleft['Y'] = upleft['Y'] - vecs.iloc[0, 1]
    upleft['Z'] = upleft['Z'] - vecs.iloc[0, 2]

    upright = up.copy()
    upright['X'] = upright['X'] + vecs.iloc[0, 0]
    upright['Y'] = upright['Y'] + vecs.iloc[0, 1]
    upright['Z'] = upright['Z'] + vecs.iloc[0, 2]

    upfront = up.copy()
    upfront['X'] = upfront['X'] - vecs.iloc[2, 0]
    upfront['Y'] = upfront['Y'] - vecs.iloc[2, 1]
    upfront['Z'] = upfront['Z'] - vecs.iloc[2, 2]

    uplfront = upfront.copy()
    uplfront['X'] = uplfront['X'] - vecs.iloc[0, 0]
    uplfront['Y'] = uplfront['Y'] - vecs.iloc[0, 1]
    uplfront['Z'] = uplfront['Z'] - vecs.iloc[0, 2]

    uprfront = upfront.copy()
    uprfront['X'] = uprfront['X'] + vecs.iloc[0, 0]
    uprfront['Y'] = uprfront['Y'] + vecs.iloc[0, 1]
    uprfront['Z'] = uprfront['Z'] + vecs.iloc[0, 2]

    upback = up.copy()
    upback['X'] = upback['X'] + vecs.iloc[2, 0]
    upback['Y'] = upback['Y'] + vecs.iloc[2, 1]
    upback['Z'] = upback['Z'] + vecs.iloc[2, 2]

    uplback = upback.copy()
    uplback['X'] = uplback['X'] - vecs.iloc[0, 0]
    uplback['Y'] = uplback['Y'] - vecs.iloc[0, 1]
    uplback['Z'] = uplback['Z'] - vecs.iloc[0, 2]

    uprback = upback.copy()
    uprback['X'] = uprback['X'] + vecs.iloc[0, 0]
    uprback['Y'] = uprback['Y'] + vecs.iloc[0, 1]
    uprback['Z'] = uprback['Z'] + vecs.iloc[0, 2]

    # Bottom row
    down = initial.copy()
    down['X'] = down['X'] - vecs.iloc[1, 0]
    down['Y'] = down['Y'] - vecs.iloc[1, 1]
    down['Z'] = down['Z'] - vecs.iloc[1, 2]

    downleft = down.copy()
    downleft['X'] = downleft['X'] - vecs.iloc[0, 0]
    downleft['Y'] = downleft['Y'] - vecs.iloc[0, 1]
    downleft['Z'] = downleft['Z'] - vecs.iloc[0, 2]

    downright = down.copy()
    downright['X'] = downright['X'] + vecs.iloc[0, 0]
    downright['Y'] = downright['Y'] + vecs.iloc[0, 1]
    downright['Z'] = downright['Z'] + vecs.iloc[0, 2]

    downfront = down.copy()
    downfront['X'] = downfront['X'] - vecs.iloc[2, 0]
    downfront['Y'] = downfront['Y'] - vecs.iloc[2, 1]
    downfront['Z'] = downfront['Z'] - vecs.iloc[2, 2]

    downlfront = downfront.copy()
    downlfront['X'] = downlfront['X'] - vecs.iloc[0, 0]
    downlfront['Y'] = downlfront['Y'] - vecs.iloc[0, 1]
    downlfront['Z'] = downlfront['Z'] - vecs.iloc[0, 2]

    downrfront = downfront.copy()
    downrfront['X'] = downrfront['X'] + vecs.iloc[0, 0]
    downrfront['Y'] = downrfront['Y'] + vecs.iloc[0, 1]
    downrfront['Z'] = downrfront['Z'] + vecs.iloc[0, 2]

    downback = down.copy()
    downback['X'] = downback['X'] + vecs.iloc[2, 0]
    downback['Y'] = downback['Y'] + vecs.iloc[2, 1]
    downback['Z'] = downback['Z'] + vecs.iloc[2, 2]

    downlback = downback.copy()
    downlback['X'] = downlback['X'] - vecs.iloc[0, 0]
    downlback['Y'] = downlback['Y'] - vecs.iloc[0, 1]
    downlback['Z'] = downlback['Z'] - vecs.iloc[0, 2]

    downrback = downback.copy()
    downrback['X'] = downrback['X'] + vecs.iloc[0, 0]
    downrback['Y'] = downrback['Y'] + vecs.iloc[0, 1]
    downrback['Z'] = downrback['Z'] + vecs.iloc[0, 2]

    super_cell = pd.concat([uplfront, upfront, uprfront, upleft, up, upright, uplback, upback, uprback,
                            lfront, front, rfront, left, initial, right, lback, back, rback,
                            downlfront, downfront, downrfront, downleft, down, downright, downlback, downback,
                            downrback],
                           axis=0, ignore_index=True)

    return super_cell


meta = pd.read_csv("jmol_xyz.mol",
                   delim_whitespace=True,
                   skiprows=5,
                   header=None,
                   nrows=1,
                   names=['A', 'B', 'C', 'Atoms', 'Bonds', 'D', 'E', 'F'])
meta = meta.drop(labels=['A', 'B', 'C', 'D', 'E', 'F'],
                 axis=1)
xyz_atoms = pd.read_csv("jmol_xyz.mol",
                        delim_whitespace=True,
                        skiprows=7,
                        header=None,
                        names=['A', 'B', 'C', 'Atom', 'X', 'Y', 'Z', 'D', 'E', 'F', 'G', 'H'],
                        nrows=meta.iloc[0, 0])
xyz_atoms = xyz_atoms.drop(labels=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                           axis=1)
xyz_bonds = pd.read_csv("jmol_xyz.mol",
                        delim_whitespace=True,
                        skiprows=7 + meta.iloc[0, 0] + 2,
                        header=None,
                        names=['A', 'B', 'C', 'D', 'Atom1', 'Atom2'],
                        nrows=meta.iloc[0, 1],
                        engine='python')
xyz_bonds = xyz_bonds.drop(labels=['A', 'B', 'C', 'D'],
                           axis=1)
vectors = pd.DataFrame([[35.4138000000, 0.0000000000, 0.0000000000],
                        [-1.7972342519, 39.7604017842, 0.0000000000],
                        [-5.7341619914, -2.8978936073, 50.9504135813]],
                       columns=['a', 'b', 'c'])
big_xyz = periodic(xyz_atoms, vectors)

big_meta = pd.read_csv("jmol_big_xyz.mol", 
                       delim_whitespace=True, skiprows=5, header=None, nrows=1,
                       names=['A', 'B', 'C', 'Atoms', 'Bonds', 'D', 'E', 'F'])
big_meta = big_meta.drop(labels=['A', 'B', 'C', 'D', 'E', 'F'], axis=1)
big_xyz_atom = pd.read_csv("jmol_big_xyz.mol", 
                           delim_whitespace=True, skiprows=7, header=None,
                           names=['A', 'B', 'C', 'Atom', 'X', 'Y', 'Z', 'D', 'E', 'F', 'G', 'H'],
                           nrows=big_meta.iloc[0, 0])
big_xyz_atom = big_xyz_atom.drop(labels=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], axis=1)
big_xyz_bonds = pd.read_csv("jmol_big_xyz.mol", 
                            delim_whitespace=True, skiprows=7 + big_meta.iloc[0, 0] + 2, header=None,
                            names=['A', 'B', 'C', 'D', 'Atom1', 'Atom2'],
                            nrows=big_meta.iloc[0, 1], engine='python')
big_xyz_bonds = big_xyz_bonds.drop(labels=['A', 'B', 'C', 'D'], axis=1)
big_xyz_bonds['Atom1'] = big_xyz_bonds['Atom1'] - 1
big_xyz_bonds['Atom2'] = big_xyz_bonds['Atom2'] - 1
bonds = big_xyz_bonds.merge(big_xyz_atom,
                            how='left',
                            left_on='Atom1',
                            right_index=True)
bonds = bonds.drop(['X', 'Y', 'Z'],
                   axis=1)
bonds = bonds.merge(big_xyz_atom,
                    how='left',
                    left_on='Atom2',
                    right_index=True)
bonds = bonds.drop(['X', 'Y', 'Z'],
                   axis=1)
bonds.columns = ['Atom1', 'Atom2', 'Atom1Name', 'Atom2Name']
G = nx.from_pandas_edgelist(big_xyz_bonds,
                            source='Atom1',
                            target='Atom2')
molecules = list(G.subgraph(c) for c in nx.connected_components(G))
molslist = []
for molecule in range(len(molecules)):
    molslist.append((molecule, list(molecules[molecule].nodes)))
slist = []
for i, n in molslist:
    molseries = pd.Series(data=i,
                          name='MOL',
                          index=n)
    slist.append(molseries)
molsdf = pd.concat(slist)
molsdf = molsdf.to_frame(name='MOL')
big_xyz_atom = big_xyz_atom.merge(molsdf,
                                  right_index=True,
                                  left_index=True,
                                  how='left')
shared = pd.merge(big_xyz_atom, xyz_atoms, how='inner', on=['X', 'Y', 'Z'])
shared.drop(['Atom_y'], axis=1, inplace=True)
shared.rename(columns={'Atom_x': 'Atom'}, inplace=True)

'''
Next....need to compare atoms that are present in both big_xyz_atom and xyz_atom. 
For molecules in big_xyz_atom, where atoms are present in xyz_atom...
need to determine the geometric centre of the molecule.

Also, should start thinking about how to set up both __init__ method / attributes, and class methods for memory 
optimisation. 

Should then think about introducing class to Soft_X-ray_Analysis.
'''
