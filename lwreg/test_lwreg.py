# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest
import sqlite3
from rdkit import Chem
from . import utils
from rdkit.Chem import RegistrationHash


class TestLWReg(unittest.TestCase):
    def setUp(self):
        if 1:
            cn = sqlite3.connect(':memory:')
            self._config = {'connection': cn}
        else:
            self._config = {'dbfile': 'test.sqlt'}

    def baseRegister(self):
        smis = ('CC[C@H](F)Cl', 'CC[C@@H](F)Cl', 'CCC(F)Cl', 'CC(F)(Cl)C')
        utils.initdb(config=self._config)
        for smi in smis:
            utils.register(smiles=smi, config=self._config)
        mols = [Chem.MolFromSmiles(x) for x in ('Cc1[nH]ncc1', 'Cc1n[nH]cc1')]
        for mol in mols:
            utils.register(mol=mol, config=self._config)

    def testRegister(self):
        utils.initdb(config=self._config)
        self.assertEqual(utils.register(smiles='CCC', config=self._config), 1)
        self.assertEqual(utils.register(smiles='CCCO', config=self._config), 2)
        self.assertRaises(
            sqlite3.IntegrityError,
            lambda: utils.register(smiles='CCC', config=self._config))
        self.assertRaises(
            sqlite3.IntegrityError,
            lambda: utils.register(smiles='CCC.O', config=self._config))
        self.assertEqual(utils.register(smiles='CCCOC', config=self._config),
                         3)
        self.assertEqual(
            utils.register(mol=Chem.MolFromSmiles('C1CCC1'),
                           config=self._config), 4)

    def testBulkRegister(self):
        utils.initdb(config=self._config)
        mols = [Chem.MolFromSmiles(x) for x in ('CCC','CCCO','C1','c1cc1','CCC','C1CC1')]
        self.assertEqual(utils.bulk_register(mols=mols, config=self._config), [1,2,None,None,None,3])

    def testQuery(self):
        self.baseRegister()
        self.assertEqual(utils.query(smiles='CC(F)Cl', config=self._config),
                         [])
        self.assertEqual(utils.query(smiles='CCC(F)Cl', config=self._config),
                         [3])
        self.assertEqual(
            utils.query(mol=Chem.MolFromSmiles('CCC(F)Cl'),
                        config=self._config), [3])
        self.assertEqual(
            utils.query(mol=Chem.MolFromSmiles('CCC(F)Cl.O'),
                        config=self._config), [3])
        self.assertEqual(
            utils.query(smiles='CCC(F)Cl',
                        layers='NO_STEREO_SMILES',
                        config=self._config), [1, 2, 3])
        self.assertEqual(
            utils.query(smiles='CCC(F)Cl',
                        layers=[RegistrationHash.HashLayer.NO_STEREO_SMILES],
                        config=self._config), [1, 2, 3])
        self.assertEqual(
            utils.query(smiles='CCC(F)Cl',
                        layers=[RegistrationHash.HashLayer.FORMULA],
                        config=self._config), [1, 2, 3, 4])

        self.assertEqual(
            utils.query(smiles='Cc1[nH]ncc1',
                        layers=[RegistrationHash.HashLayer.TAUTOMER_HASH],
                        config=self._config), [5, 6])

    def testRetrieve(self):
        self.baseRegister()
        res = utils.retrieve(ids=[12], config=self._config)
        self.assertEqual(len(res), 0)
        res = utils.retrieve(ids=[1], config=self._config)
        self.assertEqual(len(res), 1)
        tpl = res[0]
        self.assertEqual(tpl[0], 1)
        mb = Chem.MolFromMolBlock(tpl[1])
        self.assertEqual(
            utils.query(smiles=Chem.MolToSmiles(mb), config=self._config), [1])


if __name__ == '__main__':
    unittest.main()