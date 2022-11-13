# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest
import sqlite3
from rdkit import Chem
from . import utils

class TestLWReg(unittest.TestCase):
    def setUp(self):
        if 1:
            cn = sqlite3.connect(':memory:')
            self._config = {'connection':cn}
        else:
            self._config = {'dbfile':'test.sqlt'}

    def testRegister(self):
        utils.initdb(config=self._config)
        self.assertEqual(utils.register(smiles='CCC',config=self._config),1)
        self.assertEqual(utils.register(smiles='CCCO',config=self._config),2)
        self.assertRaises(sqlite3.IntegrityError,lambda : utils.register(smiles='CCC',config=self._config))
        self.assertEqual(utils.register(smiles='CCCOC',config=self._config),3)
        self.assertEqual(utils.register(mol=Chem.MolFromSmiles('C1CCC1'),config=self._config),4)




if __name__ == '__main__':
    unittest.main()