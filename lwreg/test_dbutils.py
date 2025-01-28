# Copyright (C) 2025 Greg Landrum and other lwreg contributors
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest

from rdkit import Chem

import json
import numpy as np
import copy
from collections import namedtuple

try:
    from . import utils
    from . import db_utils
except ImportError:
    import utils
    import db_utils

try:
    import psycopg2
except ImportError:
    psycopg2 = None
if psycopg2:
    # we have the connector for postgresql. Is there a server running?
    cfg = utils.defaultConfig()
    cfg['dbname'] = 'lwreg_tests'
    #cfg['host'] = 'localhost'
    cfg['dbtype'] = 'postgresql'
    try:
        cn = utils.connect(config=cfg)
        curs = cn.cursor()
        curs.execute(
            "select * from pg_available_extensions where name='rdkit'")
        if not curs.fetchone():
            rdkit_cartridge_present = False
        else:
            rdkit_cartridge_present = True
    except psycopg2.OperationalError:
        # server not running
        psycopg2 = None
        rdkit_cartridge_present = False


@unittest.skipIf(psycopg2 is None or not rdkit_cartridge_present,
                 "skipping postgresql tests")
class TestCartridgePgSQL(unittest.TestCase):
    integrityError = psycopg2.errors.UniqueViolation if psycopg2 else None

    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['dbname'] = 'lwreg_tests'
        self._config['dbtype'] = 'postgresql'

    def tearDown(self):
        utils._initdb(config=self._config, confirm=True)
        utils._clear_cached_connection()
        self._config = None
        return super().tearDown()

    def baseRegister(self):
        smis = ('CC[C@H](F)Cl', 'CC[C@@H](F)Cl', 'CCC(F)Cl', 'CC(F)(Cl)C')
        utils._initdb(config=self._config, confirm=True)
        for smi in smis:
            utils.register(smiles=smi, config=self._config)
        mols = [Chem.MolFromSmiles(x) for x in ('Cc1[nH]ncc1', 'Cc1n[nH]cc1')]
        for mol in mols:
            utils.register(mol=mol, config=self._config)

    def testPopulateSchema(self):
        self.baseRegister()
        db_utils.populate_rdkit_schema(self._config, force=True)
        conn = utils.connect(config=self._config)
        curs = conn.cursor()
        curs.execute("select count(*) from rdk.mols")
        self.assertEqual(curs.fetchone()[0], 6)
        curs = None

        # make sure the trigger is there too:
        utils.register(smiles='CCC[C@H](F)Cl', config=self._config)
        curs = conn.cursor()
        curs.execute("select count(*) from rdk.mols")
        self.assertEqual(curs.fetchone()[0], 7)
        curs = None


@unittest.skipIf(psycopg2 is None or not rdkit_cartridge_present,
                 "skipping postgresql tests")
class TestCartridgePgSQLWithSchema(TestCartridgePgSQL):

    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['dbname'] = 'lwreg_tests'
        self._config['dbtype'] = 'postgresql'
        self._config['lwregSchema'] = 'lwreg'


@unittest.skipIf(psycopg2 is None or not rdkit_cartridge_present,
                 "skipping postgresql tests")
class TestCartridgePgSQLWithRDKitSchema(TestCartridgePgSQL):

    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['dbname'] = 'lwreg_tests'
        self._config['dbtype'] = 'postgresql'
        self._config['rdkitSchema'] = 'rdkit_test'

    def testPopulateSchema(self):
        self.baseRegister()
        db_utils.populate_rdkit_schema(self._config, force=True)
        conn = utils.connect(config=self._config)
        curs = conn.cursor()
        curs.execute("select count(*) from rdkit_test.mols")
        self.assertEqual(curs.fetchone()[0], 6)
        curs = None

        # make sure the trigger is there too:
        utils.register(smiles='CCC[C@H](F)Cl', config=self._config)
        curs = conn.cursor()
        curs.execute("select count(*) from rdkit_test.mols")
        self.assertEqual(curs.fetchone()[0], 7)
        curs = None


class TestToJSON(unittest.TestCase):

    def testBasics(self):
        data = {
            'a': 1,
            'b': 'foo',
            'c': [1, 2, 3],
        }
        self.assertEqual(db_utils.dict_to_json(data), json.dumps(data))

    def testNumpy(self):
        data = {
            'a': 1,
            'b': 'foo',
            'c': np.array([1, 2, 3]),
        }
        encData = copy.deepcopy(data)
        encData['c'] = encData['c'].tolist()

        self.assertEqual(db_utils.dict_to_json(data), json.dumps(encData))

    def testNamedTuple(self):
        data = {
            'a': 1,
            'b': 'foo',
            'c': namedtuple('foo', ['a', 'b'])(1, 2),
        }
        encData = copy.deepcopy(data)
        encData['c'] = {'a': 1, 'b': 2}

        self.assertEqual(db_utils.dict_to_json(data), json.dumps(encData))

    def testKwargs(self):
        data = {
            'c': 1,
            'b': 'foo',
            'a': [1, 2, 3],
        }
        self.assertEqual(db_utils.dict_to_json(data),
                         '{"c": 1, "b": "foo", "a": [1, 2, 3]}')
        self.assertEqual(db_utils.dict_to_json(data, sort_keys=True),
                         '{"a": [1, 2, 3], "b": "foo", "c": 1}')


if __name__ == '__main__':
    unittest.main()
