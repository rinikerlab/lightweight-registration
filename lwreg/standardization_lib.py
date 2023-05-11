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


class OverlappingAtomsCheck(Standardization):
    name = "has_overlapping_atoms"
    explanation = "molecule has two (or more) atoms with coordinates within a threshold distance of each other"
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
