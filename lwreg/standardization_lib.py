# Copyright (C) 2023 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,

from rdkit import Chem


class Standardization:
    '''
    a Standardization should have a __call__ method which returns either the 
      standardized molecule (on success) or None (on failure)
    '''
    __slots__ = ["name", "explanation"]

    def __call__(self, mol):
        return mol


class OverlappingAtomsCheck(Standardization):
    name = "has_overlapping_atoms"
    explanation = "fails if molecule has at least two atoms which are closer than a threshold distance to each other"
    threshold = 0.0001

    def __call__(self, mol):
        if mol.GetNumConformers():
            t2 = self.threshold * self.threshold
            conf = mol.GetConformer()
            pts = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            for i, pti in enumerate(pts):
                for j in range(i):
                    d = pti - pts[j]
                    if d.LengthSq() < t2:
                        return None
        return mol


class PolymerCheck(Standardization):
    name = "has_polymer_info"
    explanation = "fails if molecule has an SGroup associated with polymers"
    polymerTypes = ['SRU', 'COP', 'MON', 'CRO', 'GRA']

    def __call__(self, mol):
        for sg in Chem.GetMolSubstanceGroups(mol):
            typ = sg.GetProp('TYPE')
            if typ in self.polymerTypes:
                return None

        return mol
