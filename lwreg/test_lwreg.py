# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest
import sqlite3

from . import utils

class TestLWReg(unittest.TestCase):
    def setUp(self):
        cn = sqlite3.connect(':memory:')
        self._config = {'connection':cn}

    def testRegister(self):
        utils.initdb(config=self._config)
        self.assertEqual(utils.register(smiles='CCC',config=self._config),1)
        self.assertEqual(utils.register(smiles='CCCO',config=self._config),2)
        self.assertRaises(sqlite3.IntegrityError,lambda : utils.register(smiles='CCC',config=self._config))

if __name__ == '__main__':
    unittest.main()