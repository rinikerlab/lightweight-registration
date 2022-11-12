# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE, 

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import RegistrationHash

from collections import namedtuple

MolTuple = namedtuple('MolTuple',('mol','datatype','rawdata'))
def parse_mol(molfile=None, molblock=None, smiles=None, config=None):
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        datatype = 'smiles'
        raw = smiles
    return MolTuple(mol,datatype,raw)
    
def standardize_mol(mol, config=None):
    sMol = rdMolStandardize.FragmentParent(mol)
    return sMol

def hash_mol(mol,escape=None,config=None):
    layers = RegistrationHash.GetMolLayers(mol,escape=escape)
    mhash = RegistrationHash.GetMolHash(layers)
    return mhash,layers
