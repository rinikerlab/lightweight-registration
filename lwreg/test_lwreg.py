# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import os
import time
from datetime import datetime, timedelta
import unittest
import sqlite3
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolTransforms
import random
import copy
import tempfile

try:
    from . import utils
    from .utils import RegistrationFailureReasons
    from . import standardization_lib
except ImportError:
    import utils
    from utils import RegistrationFailureReasons
    import standardization_lib

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
    except psycopg2.OperationalError:
        # server not running
        psycopg2 = None


class TestLWReg(unittest.TestCase):
    integrityError = sqlite3.IntegrityError

    def setUp(self):
        cn = sqlite3.connect(':memory:')
        self._config = utils.defaultConfig()
        self._config['connection'] = cn

    def baseRegister(self):
        smis = ('CC[C@H](F)Cl', 'CC[C@@H](F)Cl', 'CCC(F)Cl', 'CC(F)(Cl)C')
        utils._initdb(config=self._config, confirm=True)
        for smi in smis:
            utils.register(smiles=smi, config=self._config)
        mols = [Chem.MolFromSmiles(x) for x in ('Cc1[nH]ncc1', 'Cc1n[nH]cc1')]
        for mol in mols:
            utils.register(mol=mol, config=self._config)

    def testRegister(self):
        utils._initdb(config=self._config, confirm=True)

        self.assertEqual(utils.registration_counts(config=self._config), 0)
        self.assertEqual(utils.get_all_registry_numbers(config=self._config),
                         ())

        self.assertEqual(utils.register(smiles='CCC', config=self._config), 1)
        self.assertEqual(utils.register(smiles='CCCO', config=self._config), 2)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCC', config=self._config))
        self.assertEqual(
            utils.register(smiles='CCCO',
                           config=self._config,
                           fail_on_duplicate=False), 2)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCC.O', config=self._config))
        self.assertGreater(utils.register(smiles='CCCOC', config=self._config),
                           2)
        self.assertGreater(
            utils.register(mol=Chem.MolFromSmiles('C1CCC1'),
                           config=self._config), 3)
        self.assertEqual(utils.register(mol=None, config=self._config),
                         RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(utils.register(smiles='CCCC1', config=self._config),
                         RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(utils.register(smiles='c1nccc1', config=self._config),
                         RegistrationFailureReasons.FILTERED)

        self.assertEqual(utils.registration_counts(config=self._config), 4)
        expected = {
            'sqlite3': (1, 2, 3, 4),
            'postgresql': (1, 2, 6, 7),
        }
        self.assertEqual(utils.get_all_registry_numbers(config=self._config),
                         expected[self._config['dbtype']])

    def testGetDelimiter(self):
        self.assertEqual(
            utils._get_delimiter('test_data/test_smiles_no_delim.smi'), None)
        self.assertEqual(
            utils._get_delimiter(
                'test_data/test_smiles_no_delim_with_header.smi'), None)
        self.assertEqual(
            utils._get_delimiter('test_data/test_smiles_with_header.smi'), ';')
        self.assertEqual(utils._get_delimiter('test_data/test_smiles.smi'),
                         ' ')

    def testGetSmilesColumn(self):
        self.assertEqual(
            utils._get_smiles_column('test_data/test_smiles_with_header.smi',
                                     delimiter=';'), 1)
        self.assertEqual(
            utils._get_smiles_column('test_data/test_smiles.smi',
                                     delimiter=' '), 1)

    def testHasHeader(self):
        self.assertTrue(
            utils._has_header('test_data/test_smiles_no_delim_with_header.smi',
                              delimiter=' '))
        self.assertFalse(
            utils._has_header('test_data/test_smiles_no_delim.smi',
                              delimiter=' '))
        self.assertTrue(
            utils._has_header('test_data/test_smiles_with_header.smi',
                              delimiter=','))
        self.assertFalse(
            utils._has_header('test_data/test_smiles.smi', delimiter=' '))

    def testGetMolsFromSmilesfile(self):
        filenames = [
            'test_data/test_smiles_no_delim.smi',
            'test_data/test_smiles_no_delim_with_header.smi',
            'test_data/test_smiles_with_header.smi',
            'test_data/test_smiles.smi'
        ]
        for filename in filenames:
            mols = utils._get_mols_from_smilesfile(filename)
            self.assertEqual(len(mols), 6)

    def testBulkRegister(self):
        utils._initdb(config=self._config, confirm=True)
        mols = [
            Chem.MolFromSmiles(x)
            for x in ('CCC', 'CCCO', 'C1', 'c1cc1', 'CCC', 'C1CC1')
        ]
        res = utils.bulk_register(mols=mols, config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res[2], RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(res[3], RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(res[4], RegistrationFailureReasons.DUPLICATE)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         2)
        self.assertEqual(res.count(RegistrationFailureReasons.DUPLICATE), 1)

        res = utils.bulk_register(sdfile='test_data/test_molecules.sdf',
                                  config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         0)
        self.assertEqual(res.count(RegistrationFailureReasons.DUPLICATE), 0)

        res = utils.bulk_register(
            smilesfile='test_data/test_smiles_no_delim_with_header.smi',
            config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         0)
        self.assertEqual(res.count(RegistrationFailureReasons.DUPLICATE), 0)

        res = utils.bulk_register(
            smilesfile='test_data/test_smiles_with_header.smi',
            config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         0)
        self.assertEqual(res.count(RegistrationFailureReasons.DUPLICATE), 0)

        res = utils.bulk_register(
            smilesfile='test_data/test_smiles_no_delim.smi',
            config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         0)
        self.assertEqual(res.count(RegistrationFailureReasons.DUPLICATE), 0)

        res = utils.bulk_register(smilesfile='test_data/test_smiles.smi',
                                  config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         0)
        self.assertEqual(res.count(RegistrationFailureReasons.DUPLICATE), 0)
        ## add test for config where Hs are to be kept and molecules come from SD file
        configWithHs = self._config
        configWithHs["removeHs"] = 0
        configWithHs["standardization"] = 'none'
        utils._initdb(config=configWithHs, confirm=True)
        filename = 'test_data/test_molecules_hs.sdf'
        res = utils.bulk_register(sdfile=filename,
                                  config=configWithHs)
        res = utils.retrieve(ids=[1], config=configWithHs)
        smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(res[1][0],removeHs=False))
        self.assertEqual(smiles,'[H]C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]')

    def testBulkRegisterAllowDupes(self):
        utils._initdb(config=self._config, confirm=True)
        mols = [
            Chem.MolFromSmiles(x)
            for x in ('CCC', 'CCCO', 'C1', 'c1cc1', 'CCC', 'C1CC1')
        ]

        res = utils.bulk_register(mols=mols,
                                  fail_on_duplicate=False,
                                  config=self._config)
        self.assertEqual(len(res), 6)
        self.assertEqual(res[2], RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(res[3], RegistrationFailureReasons.PARSE_FAILURE)
        self.assertNotIn(RegistrationFailureReasons.DUPLICATE, res)
        self.assertEqual(res.count(RegistrationFailureReasons.PARSE_FAILURE),
                         2)

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
                        layers=[utils.HashLayer.NO_STEREO_SMILES],
                        config=self._config), [1, 2, 3])
        self.assertEqual(
            utils.query(smiles='CCC(F)Cl',
                        layers=utils.HashLayer.NO_STEREO_SMILES,
                        config=self._config), [1, 2, 3])
        self.assertEqual(
            utils.query(smiles='CCC(F)Cl',
                        layers=[utils.HashLayer.FORMULA],
                        config=self._config), [1, 2, 3, 4])

        self.assertEqual(
            utils.query(smiles='Cc1[nH]ncc1',
                        layers=[utils.HashLayer.TAUTOMER_HASH],
                        config=self._config), [5, 6])

    def testRetrieve(self):
        self.baseRegister()
        res = utils.retrieve(ids=[12], config=self._config)
        self.assertEqual(len(res), 0)
        res = utils.retrieve(ids=[1], config=self._config)
        self.assertEqual(len(res), 1)
        tpl = res[1]
        self.assertEqual(len(tpl), 2)
        mb = Chem.MolFromMolBlock(tpl[0])
        self.assertEqual(
            utils.query(smiles=Chem.MolToSmiles(mb), config=self._config), [1])
        self.assertEqual(tpl[1], 'mol')
        res = utils.retrieve(ids=[100], config=self._config)
        self.assertEqual(len(res), 0)

        res = utils.retrieve(ids=[1, 2], config=self._config, as_hashes=True)
        for k, v in res.items():
            self.assertFalse('molregno' in v)
            self.assertTrue('fullhash' in v)
            self.assertTrue('canonical_smiles' in v)

        res = utils.retrieve(ids=[1, 5],
                             config=self._config,
                             as_submitted=True)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[1][1], 'smiles')
        self.assertEqual(res[5][1], 'pkl')
        m = Chem.Mol(res[5][0])
        self.assertEqual(Chem.MolToSmiles(m), 'Cc1ccn[nH]1')

    def testPickleBug(self):
        ipkl = r'\xefbeadde00000000100000000100000000000000030000000200000080010600600000000103060060000000020208006000000001010b00010001020042000000001709000000000000003f00000000126400000003000f0000005f5f636f6d707574656450726f7073060200000000000000070000006e756d41726f6d0f0000005f53746572656f6368656d446f6e65070000006e756d41726f6d01000000000f0000005f53746572656f6368656d446f6e650101000000133ab400000002000f0000005f5f636f6d707574656450726f7073060100000000000000080000005f43495052616e6b080000005f43495052616e6b02000000000002000f0000005f5f636f6d707574656450726f7073060100000000000000080000005f43495052616e6b080000005f43495052616e6b02010000000002000f0000005f5f636f6d707574656450726f7073060100000000000000080000005f43495052616e6b080000005f43495052616e6b0202000000001316'
        npkl = utils._parsePickleFromDB(ipkl)
        m = Chem.Mol(npkl)
        self.assertEqual(Chem.MolToSmiles(m), 'CCO')

    def testStandardizationOptions(self):
        lconfig = self._config.copy()
        lconfig['standardization'] = 'charge'
        utils._initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='CCCO.[Na]', config=lconfig), 1)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCC[O-]', config=lconfig))

        lconfig['standardization'] = 'tautomer'
        utils._initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='Cc1[nH]ncc1', config=lconfig),
                         1)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='Cc1n[nH]cc1', config=lconfig))

        lconfig['standardization'] = 'ChargeParent'
        utils._initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='CCCO.[Na]', config=lconfig), 1)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCC[O-]', config=lconfig))

    def testStandardizationFunctions(self):
        lconfig = self._config.copy()

        def evenAtomCount(mol):
            if mol.GetNumAtoms() % 2:
                return None
            return mol

        # Silly example which only accepts even numbers of atoms
        lconfig['standardization'] = evenAtomCount

        utils._initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='CCCO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='CCC', config=lconfig),
                         RegistrationFailureReasons.FILTERED)

        utils._initdb(config=lconfig, confirm=True)
        checker = standardization_lib.OverlappingAtomsCheck()
        lconfig['standardization'] = checker
        self.assertEqual(
            utils.register(smiles='CC |(0,0,;0,1,)|', config=lconfig), 1)
        self.assertEqual(
            utils.register(smiles='CO |(0,0,;0,0,)|', config=lconfig),
            RegistrationFailureReasons.FILTERED)

        def hasAnOxygen(mol):
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 8:
                    return mol
            return None

        # ----------------------
        # chaining standardization functions
        utils._initdb(config=lconfig, confirm=True)
        lconfig['standardization'] = [evenAtomCount, hasAnOxygen]
        self.assertEqual(utils.register(smiles='CC', config=lconfig),
                         RegistrationFailureReasons.FILTERED)
        self.assertEqual(utils.register(smiles='CO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='CCO', config=lconfig),
                         RegistrationFailureReasons.FILTERED)

        # ----------------------
        # mixing standardization options and functions
        utils._initdb(config=lconfig, confirm=True)
        lconfig['standardization'] = ('fragment', evenAtomCount, hasAnOxygen)
        self.assertEqual(utils.register(smiles='CC', config=lconfig),
                         RegistrationFailureReasons.FILTERED)
        self.assertEqual(utils.register(smiles='CO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='C[O-].[Na]', config=lconfig),
                         2)
        self.assertEqual(utils.register(smiles='CCO', config=lconfig),
                         RegistrationFailureReasons.FILTERED)

    def testStandardizationLibCheckers(self):
        lconfig = self._config.copy()

        utils._initdb(config=lconfig, confirm=True)
        checker = standardization_lib.OverlappingAtomsCheck()
        lconfig['standardization'] = checker
        self.assertEqual(
            utils.register(smiles='CC |(0,0,;0,1,)|', config=lconfig), 1)
        self.assertEqual(
            utils.register(smiles='CO |(0,0,;0,0,)|', config=lconfig),
            RegistrationFailureReasons.FILTERED)

        checker = standardization_lib.PolymerCheck()
        lconfig['standardization'] = checker
        self.assertEqual(
            utils.register(smiles='*-CO-* |$star_e;;;star_e$,Sg:n:1,2::ht|',
                           config=lconfig),
            RegistrationFailureReasons.FILTERED)
        self.assertEqual(
            utils.register(molblock='''
  Mrv2305 05112307022D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.6667 0.7083 0 0
M  V30 2 C -5.333 1.4783 0 0
M  V30 3 O -3.9993 0.7083 0 0
M  V30 4 C -2.6656 1.4783 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(2 2 3) SAP=(3 2 1 1) SAP=(3 3 4 2) XBONDS=(2 1 3) -
M  V30 LABEL=CO ESTATE=E
M  V30 END SGROUP
M  V30 END CTAB
M  END
''',
                           config=lconfig), 2)

    def testTimestamp(self):
        utils._initdb(config=self._config, confirm=True)
        utils.register(smiles='CCC', config=self._config)
        time.sleep(2)
        utils.register(smiles='CCCC', config=self._config)
        cn = utils.connect(config=self._config)
        curs = cn.cursor()
        curs.execute(
            f"select molregno, timestamp from {utils.origDataTableName} order by molregno asc"
        )
        d = curs.fetchall()
        curs = None
        timestamps = []
        for row in d:
            timestamps.append(datetime.strptime(row[1], "%Y-%m-%d %H:%M:%S"))
        self.assertEqual(timestamps[1] - timestamps[0] > timedelta(0), True)

    def testConfigFromDatabase(self):
        lconfig = self._config.copy()
        lconfig['standardization'] = 'charge'
        utils._initdb(config=lconfig, confirm=True)
        self.assertEqual(
            utils.register(smiles='CCC[O-].[Na+]', config=lconfig), 1)
        self.assertEqual(
            utils.register(smiles='CCC(=O)[O-].[Na+]', config=lconfig), 2)
        self.assertEqual(utils.registration_counts(config=lconfig), 2)
        conn = None
        if 'connection' in lconfig:
            conn = lconfig['connection']
        nconfig = utils.configure_from_database(
            connection=conn,
            dbname=lconfig['dbname'],
            dbtype=lconfig['dbtype'],
            lwregSchema=lconfig['lwregSchema'])
        configCopy = lconfig.copy()
        nconfigCopy = nconfig.copy()
        if 'connection' in configCopy:
            del configCopy['connection']
            del nconfigCopy['connection']
        self.assertEqual(nconfigCopy, configCopy)

        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCC[O-]', config=nconfig))
        self.assertGreater(
            utils.register(smiles='CCCC(=O)[O-].[Na+]', config=nconfig), 2)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCCC(=O)O', config=nconfig))

    def testConfigFromDatabaseWithoutDbType(self):
        # this test only makes sense for sqlite
        if self._config['dbtype'] == 'postgresql':
            return
        tmpfile = tempfile.NamedTemporaryFile()
        tmpfile.close()
        lconfig = self._config.copy()
        if 'connection' in lconfig:
            del lconfig['connection']
        lconfig['dbname'] = tmpfile.name
        lconfig['standardization'] = 'charge'
        utils._initdb(config=lconfig, confirm=True)
        self.assertEqual(
            utils.register(smiles='CCC[O-].[Na+]', config=lconfig), 1)
        self.assertEqual(
            utils.register(smiles='CCC(=O)[O-].[Na+]', config=lconfig), 2)
        self.assertEqual(utils.registration_counts(config=lconfig), 2)
        nconfig = utils.configure_from_database(
            connection=None,
            dbname=lconfig['dbname'],
            dbtype=None,
            lwregSchema=lconfig['lwregSchema'])
        print(lconfig)
        print(nconfig)
        self.assertEqual(nconfig['dbtype'], 'sqlite3')
        self.assertEqual(nconfig, lconfig)

    def testSetDefaultConfig(self):
        lconfig = self._config.copy()
        lconfig['standardization'] = 'charge'
        utils.set_default_config(lconfig)
        utils._initdb(confirm=True)
        self.assertEqual(utils.register(smiles='CCC[O-].[Na+]'), 1)
        self.assertEqual(utils.register(smiles='CCC(=O)[O-].[Na+]'), 2)
        self.assertEqual(utils.registration_counts(), 2)
        self.assertRaises(self.integrityError,
                          lambda: utils.register(smiles='CCC(=O)O'))

    def testDbIntegrityConstraints(self):
        cfg = self._config
        utils.set_default_config(cfg)
        utils._initdb(confirm=True)
        self.assertEqual(utils.register(smiles='CCO'), 1)
        self.assertEqual(utils.register(smiles='CCOC'), 2)
        cn = utils.connect(cfg)
        curs = cn.cursor()
        curs.execute('PRAGMA foreign_keys=ON')
        self.assertRaises(
            sqlite3.IntegrityError, lambda: curs.execute(
                "insert into orig_data (molregno, data, datatype) values (4, 'foo', 'bar')"
            ))
        cn.rollback()
        self.assertRaises(
            sqlite3.IntegrityError, lambda: curs.execute(
                "insert into molblocks values (4, 'foo', 'bar')"))
        cn.rollback()


class TestLWRegTautomerv2(unittest.TestCase):
    integrityError = sqlite3.IntegrityError

    def setUp(self):
        cn = sqlite3.connect(':memory:')
        self._config = utils.defaultConfig()
        self._config['connection'] = cn
        self._config['useTautomerHashv2'] = 1

    def testRegister(self):
        utils._initdb(config=self._config, confirm=True)
        self.assertEqual(utils.register(smiles='CCC', config=self._config), 1)
        self.assertEqual(utils.register(smiles='CCCO', config=self._config), 2)
        self.assertEqual(utils.register(smiles='CC=CO', config=self._config),
                         3)
        self.assertEqual(utils.register(smiles='CCC=O', config=self._config),
                         4)

        self.assertEqual(
            utils.query(smiles='CC=CO',
                        layers=[utils.HashLayer.TAUTOMER_HASH],
                        config=self._config), [3, 4])
        self.assertEqual(
            utils.query(smiles='CC=CO',
                        layers=[utils.HashLayer.CANONICAL_SMILES],
                        config=self._config), [3])


@unittest.skipIf(psycopg2 is None, "skipping postgresql tests")
class TestLWRegPSQL(TestLWReg):
    integrityError = psycopg2.errors.UniqueViolation if psycopg2 else None

    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['dbname'] = 'lwreg_tests'
        self._config['dbtype'] = 'postgresql'

    def testTimestamp(self):
        utils._initdb(config=self._config, confirm=True)
        utils.register(smiles='CCC', config=self._config)
        time.sleep(2)
        utils.register(smiles='CCCC', config=self._config)
        cn = utils.connect(config=self._config)
        curs = cn.cursor()
        curs.execute(
            f"select molregno,timestamp from {utils.origDataTableName} order by molregno asc"
        )
        d = curs.fetchall()
        curs = None
        timestamps = []
        for row in d:
            timestamps.append(row[1])
        self.assertEqual(timestamps[1] - timestamps[0] > timedelta(0), True)

    def testDbIntegrityConstraints(self):
        cfg = self._config
        utils.set_default_config(cfg)
        utils._initdb(confirm=True)
        self.assertEqual(utils.register(smiles='CCO'), 1)
        self.assertEqual(utils.register(smiles='CCOC'), 2)
        cn = utils.connect(cfg)
        curs = cn.cursor()

        self.assertRaises(
            psycopg2.errors.ForeignKeyViolation, lambda: curs.execute(
                "insert into orig_data (molregno, data, datatype) values (4, 'foo', 'bar')"
            ))
        cn.rollback()

        self.assertRaises(
            psycopg2.errors.ForeignKeyViolation, lambda: curs.execute(
                "insert into molblocks values (4, 'foo', 'bar')"))
        cn.rollback()


@unittest.skipIf(psycopg2 is None, "skipping postgresql tests")
class TestLWRegPSQLWithSchema(TestLWRegPSQL):

    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['dbname'] = 'lwreg_tests'
        self._config['dbtype'] = 'postgresql'
        self._config['lwregSchema'] = 'lwreg'

    def testSchema(self):
        utils._initdb(config=self._config, confirm=True)
        self.assertEqual(utils.registrationMetadataTableName,
                         'lwreg.registration_metadata')


class TestStandardizationLabels(unittest.TestCase):

    def setUp(self):
        cn = sqlite3.connect(':memory:')
        self._config = utils.defaultConfig()
        self._config['connection'] = cn

    def testStandards(self):
        cfg = self._config
        for k in utils.standardizationOptions:
            cfg['standardization'] = k
            lbl = utils._get_standardization_label(cfg)
            self.assertEqual(lbl, k)

    def testCombined(self):
        cfg = self._config
        for k in utils.standardizationOptions:
            cl = ['foo', k]
            cfg['standardization'] = cl
            lbl = utils._get_standardization_label(cfg)
            self.assertEqual(lbl, '|'.join(cl))

    def testOthers(self):

        def func1(x):
            pass

        cfg = self._config
        for k in utils.standardizationOptions:
            for func in (func1, lambda x: x):
                cl = [k, func]
                cfg['standardization'] = cl
                lbl = utils._get_standardization_label(cfg)
                self.assertEqual(lbl, f'{k}|unknown')
            for func in (
                    standardization_lib.OverlappingAtomsCheck,
                    standardization_lib.PolymerCheck,
            ):
                cl = [k, func]
                cfg['standardization'] = cl
                lbl = utils._get_standardization_label(cfg)
                self.assertEqual(lbl, f'{k}|{func.name}')

    def testRecording(self):
        cfg = self._config
        utils._initdb(config=cfg, confirm=True)
        self.assertEqual(utils.register(smiles='CCO', config=cfg), 1)
        self.assertEqual(utils.register(smiles='CCOC', config=cfg), 2)
        oac = standardization_lib.OverlappingAtomsCheck()
        cfg['standardization'] = ['fragment', oac]
        self.assertEqual(utils.register(smiles='CCN', config=cfg), 3)
        self.assertEqual(utils.register(smiles='CCNC', config=cfg), 4)
        cn = utils.connect(cfg)
        curs = cn.cursor()
        curs.execute(
            f"select count(*) from {utils.molblocksTableName} where standardization is not null"
        )
        self.assertEqual(curs.fetchone()[0], 2)
        curs.execute(
            f"select count(*) from {utils.molblocksTableName} where standardization is null"
        )
        self.assertEqual(curs.fetchone()[0], 2)

    def testStandardizeMolFunction(self):
        cfg = self._config
        cfg['standardization'] = 'charge'
        utils._initdb(config=cfg, confirm=True)
        m = Chem.MolFromSmiles('CC[O-].[Na+]')
        nm = utils.standardize_mol(m, config=cfg)
        self.assertEqual(Chem.MolToSmiles(nm), 'CCO')


class TestRegisterConformers(unittest.TestCase):
    integrityError = sqlite3.IntegrityError

    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['registerConformers'] = True
        self._mol1 = Chem.AddHs(Chem.MolFromSmiles('OC(=O)CCCC'))
        rdDistGeom.EmbedMolecule(self._mol1, randomSeed=0xf00d)
        self._mol2 = Chem.Mol(self._mol1)
        rdDistGeom.EmbedMolecule(self._mol2, randomSeed=0xf00d + 1)
        self._mol3 = Chem.AddHs(Chem.MolFromSmiles('CCOC(=O)CCCC'))
        rdDistGeom.EmbedMolecule(self._mol3, randomSeed=0xf00d)

        m1 = Chem.MolFromSmiles(
            'F[C](Cl)(Br)I |(-0.212215,0.280139,1.44255;-0.0440094,0.0353471,0.1035;-0.457136,1.51641,-0.774943;-1.25605,-1.35502,-0.457339;1.96941,-0.476872,-0.313764)|'
        )
        m2 = Chem.MolFromSmiles(
            'F[C](Cl)(Br)I |(-0.265874,-0.363334,1.43723;-0.050665,-0.040347,0.110315;-0.530745,-1.4506,-0.855487;-1.15413,1.45519,-0.417978;2.00142,0.399096,-0.274081)|'
        )
        self.assertNotEqual(Chem.MolToSmiles(m1), Chem.MolToSmiles(m2))
        m1.AddConformer(m2.GetConformer(), assignId=True)
        self._chiralMol = m1

    def testConformerDupes(self):
        utils._initdb(config=self._config, confirm=True)
        self.assertEqual(utils.register(mol=self._mol1, config=self._config),
                         (1, 1))
        self.assertEqual(utils.register(mol=self._mol2, config=self._config),
                         (1, 2))

        aorder = list(range(self._mol1.GetNumAtoms()))
        random.shuffle(aorder)
        nmol = Chem.RenumberAtoms(self._mol1, aorder)

        # make sure we fail if it's a conformer dupe:
        with self.assertRaises(self.integrityError):
            utils.register(mol=nmol, config=self._config)

        # conformer dupe with dupes allowed:
        self.assertEqual(
            utils.register(mol=nmol,
                           config=self._config,
                           fail_on_duplicate=False), (1, 1))

    def testBulkConformers(self):
        utils._initdb(config=self._config, confirm=True)
        aorder = list(range(self._mol1.GetNumAtoms()))
        random.shuffle(aorder)
        nmol = Chem.RenumberAtoms(self._mol1, aorder)
        expected = {
            'sqlite3': ((1, 1), (1, 2), (1, 1), (2, 3)),
            'postgresql': ((1, 1), (1, 2), (1, 1), (4, 4)),
        }
        self.assertEqual(
            utils.bulk_register(mols=(self._mol1, self._mol2, nmol,
                                      self._mol3),
                                fail_on_duplicate=False,
                                config=self._config),
            expected[self._config['dbtype']])
        self.assertEqual(utils.registration_counts(config=self._config),
                         (2, 3))
        expected = {
            'sqlite3': (1, 2),
            'postgresql': (1, 4),
        }
        self.assertEqual(utils.get_all_registry_numbers(config=self._config),
                         expected[self._config['dbtype']])

    def testNoConformers(self):
        utils._initdb(config=self._config, confirm=True)
        with self.assertRaises(ValueError):
            utils.register(smiles='c1ccccc1', config=self._config)

        # we can register "empty" conformers:
        m = Chem.MolFromSmiles('c1ccccc1')
        conf = Chem.Conformer(m.GetNumAtoms())
        m.AddConformer(conf)
        self.assertEqual(utils.register(mol=m, config=self._config), (1, 1))

        # but of course the second time is a duplicate
        m = Chem.MolFromSmiles('c1ccccc1')
        conf = Chem.Conformer(m.GetNumAtoms())
        m.AddConformer(conf)
        with self.assertRaises(self.integrityError):
            utils.register(mol=m, config=self._config)

    def testMultiConfMolecule(self):
        utils._initdb(config=self._config, confirm=True)

        mol = Chem.Mol(self._mol1)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 10, randomSeed=0xf00d)
        self.assertEqual(len(cids), 10)
        # add a duplicate conformer to ensure that is handled correctly
        mol.AddConformer(mol.GetConformer())

        rres = utils.register_multiple_conformers(mol=mol,
                                                  fail_on_duplicate=False,
                                                  config=self._config)
        self.assertEqual(len(rres), 11)
        self.assertEqual(len(set(rres)), 10)
        self.assertEqual(len(set([mrn for mrn, cid in rres])), 1)

        # make sure we can add more conformers:
        mol2 = Chem.Mol(mol)
        cids = rdDistGeom.EmbedMultipleConfs(mol2, 10, randomSeed=0xd00f)
        self.assertEqual(len(cids), 10)
        rres = utils.register_multiple_conformers(mol=mol2,
                                                  fail_on_duplicate=True,
                                                  config=self._config)
        self.assertEqual(len(rres), 10)
        self.assertEqual(len(set(rres)), 10)
        self.assertEqual(len(set([mrn for mrn, cid in rres])), 1)

        # make sure we can still fail on duplicate conformers:
        utils._initdb(config=self._config, confirm=True)
        with self.assertRaises(self.integrityError):
            utils.register_multiple_conformers(mol=mol,
                                               fail_on_duplicate=True,
                                               config=self._config)

    def testConformerQuery(self):
        ''' querying using a molecule which has conformers '''
        utils._initdb(config=self._config, confirm=True)
        regids = utils.bulk_register(mols=(self._mol1, self._mol3),
                                     config=self._config)
        self.assertEqual(
            sorted(utils.query(mol=self._mol1, config=self._config)),
            [regids[0]])
        # matches topology, but not conformer
        self.assertEqual(
            sorted(utils.query(mol=self._mol2, config=self._config)), [])
        self.assertEqual(
            sorted(utils.query(mol=self._mol3, config=self._config)),
            [regids[1]])

        # query with no conformer
        qm = Chem.Mol(self._mol1)
        qm.RemoveAllConformers()
        self.assertEqual(sorted(utils.query(mol=qm, config=self._config)),
                         [regids[0][0]])

    def testConformerRetrieve(self):
        ''' querying using a molecule which has conformers '''
        utils._initdb(config=self._config, confirm=True)
        regids = utils.bulk_register(mols=(self._mol1, self._mol2, self._mol3),
                                     config=self._config)

        res = utils.retrieve(ids=(regids[0], regids[2]), config=self._config)
        self.assertEqual(len(res), 2)
        self.assertTrue(regids[0] in res)
        self.assertTrue(regids[2] in res)
        self.assertIn('M  END', res[regids[0]][0])
        self.assertEqual(res[regids[0]][1], 'mol')
        self.assertIn('M  END', res[regids[2]][0])
        self.assertEqual(res[regids[2]][1], 'mol')

        # query with just molregnos... then we get back the same thing as if we
        # were not in conformer mode.
        mrn0 = regids[0][0]
        mrn2 = regids[2][0]

        res = utils.retrieve(id=mrn0, config=self._config)
        self.assertIn(mrn0, res)
        self.assertIn('M  END', res[mrn0][0])

        res = utils.retrieve(ids=(mrn0, mrn2), config=self._config)
        self.assertEqual(len(res), 2)
        self.assertTrue(mrn0 in res)
        self.assertTrue(mrn2 in res)
        self.assertIn('M  END', res[mrn0][0])
        self.assertEqual(res[mrn0][1], 'mol')
        self.assertIn('M  END', res[mrn2][0])
        self.assertEqual(res[mrn2][1], 'mol')

        res = utils.retrieve(ids=(mrn0, mrn2),
                             config=self._config,
                             as_hashes=True)
        self.assertIn(mrn0, res)
        self.assertIn(mrn2, res)
        self.assertIn('fullhash', res[mrn0])

    def testConformerQueryById(self):
        utils._initdb(config=self._config, confirm=True)
        regids = utils.bulk_register(mols=(self._mol1, self._mol2, self._mol3),
                                     config=self._config)
        mrns, cids = zip(*regids)
        self.assertEqual(
            sorted(utils.query(ids=mrns[0:1], config=self._config)), [(1, 1),
                                                                      (1, 2)])
        expected = {
            'sqlite3': [(1, 1), (1, 2), (2, 3)],
            'postgresql': [(1, 1), (1, 2), (3, 3)],
        }
        self.assertEqual(sorted(utils.query(ids=mrns, config=self._config)),
                         expected[self._config['dbtype']])
        self.assertEqual(
            sorted(utils.query(ids=tuple(reversed(mrns)),
                               config=self._config)),
            expected[self._config['dbtype']])
        with self.assertRaises(ValueError):
            cnf = copy.deepcopy(self._config)
            cnf['registerConformers'] = False
            utils.query(ids=mrns, config=cnf)

    def testBulkConformersAndChirality(self):
        utils._initdb(config=self._config, confirm=True)

        expected = {
            'sqlite3': (
                (1, 1),
                (2, 2),
            ),
            'postgresql': (
                (1, 1),
                (2, 2),
            ),
        }
        self.assertEqual(
            utils.register_multiple_conformers(mol=self._chiralMol,
                                               fail_on_duplicate=False,
                                               config=self._config),
            expected[self._config['dbtype']])
        self.assertEqual(utils.registration_counts(config=self._config),
                         (2, 2))

    def testRegisterMolWithConfId(self):
        utils._initdb(config=self._config, confirm=True)
        self.assertEqual(
            utils.register(mol=self._chiralMol,
                           fail_on_duplicate=False,
                           config=self._config), (1, 1))
        self.assertEqual(
            utils.register(mol=self._chiralMol,
                           fail_on_duplicate=False,
                           confId=1,
                           config=self._config), (2, 2))

        self.assertEqual(
            utils.register(mol=self._chiralMol,
                           fail_on_duplicate=False,
                           confId=0,
                           config=self._config), (1, 1))

    def testRegisterThenBulkRegister(self):
        utils._initdb(config=self._config, confirm=True)
        self.assertEqual(
            utils.register(mol=self._chiralMol,
                           fail_on_duplicate=False,
                           confId=1,
                           config=self._config), (1, 1))

        expected = {
            'sqlite3': (
                (2, 2),
                (1, 1),
            ),
            'postgresql': (
                (2, 2),
                (1, 1),
            ),
        }
        self.assertEqual(
            utils.register_multiple_conformers(mol=self._chiralMol,
                                               fail_on_duplicate=False,
                                               config=self._config),
            expected[self._config['dbtype']])

    def testRegisterThenBulkRegister2(self):
        utils._initdb(config=self._config, confirm=True)
        self.assertEqual(
            utils.register(mol=self._chiralMol,
                           fail_on_duplicate=False,
                           config=self._config), (1, 1))

        expected = {
            'sqlite3': (
                (1, 1),
                (2, 2),
            ),
            'postgresql': (
                (1, 1),
                (3, 3),
            ),
        }
        self.assertEqual(
            utils.register_multiple_conformers(mol=self._chiralMol,
                                               fail_on_duplicate=False,
                                               config=self._config),
            expected[self._config['dbtype']])
        
    def testConformerStandardization(self):
        cfg = self._config.copy()
        cfg['standardization'] = [
            standardization_lib.CanonicalizeOrientation()
        ]
        utils._initdb(config=cfg, confirm=True)
        self.assertEqual(utils.register(mol=self._mol1, config=cfg), (1, 1))
        cp = Chem.Mol(self._mol1)
        conf = cp.GetConformer()
        for i in range(conf.GetNumAtoms()):
            pi = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pi.x + 0.5, pi.y - 0.3, pi.z + 1.5))
        self.assertRaises(self.integrityError,
                          lambda: utils.register(mol=cp, config=cfg))

    def testMultiConformerStandardization(self):
        cfg = self._config.copy()
        cfg['standardization'] = [
            standardization_lib.CanonicalizeOrientation()
        ]
        utils._initdb(config=cfg, confirm=True)
        cp = Chem.Mol(self._mol1)
        conf = cp.GetConformer()
        for i in range(conf.GetNumAtoms()):
            pi = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pi.x + 0.5, pi.y - 0.3, pi.z + 1.5))
        cp.AddConformer(self._mol1.GetConformer(), assignId=True)
        self.assertRaises(
            self.integrityError, lambda: utils.register_multiple_conformers(
                mol=cp, fail_on_duplicate=True, config=cfg))


@unittest.skipIf(psycopg2 is None, "skipping postgresql tests")
class TestRegisterConformersPSQL(TestRegisterConformers):
    integrityError = psycopg2.errors.UniqueViolation if psycopg2 else None

    def setUp(self):
        super(TestRegisterConformersPSQL, self).setUp()
        #getlogin = lambda: pwd.getpwuid(os.getuid())[0]
        self._config['dbname'] = 'lwreg_tests'
        self._config['dbtype'] = 'postgresql'
        self._config['password'] = 'testpw'
        #self._config['user'] = getlogin()

    def testNoSecretsInRegistrationMetadata(self):
        """Make sure initdb is not storing any secrets."""
        utils._initdb(config=self._config, confirm=True)
        with utils.connect(self._config).cursor() as cursor:
            cursor.execute('select * from registration_metadata;')
            keys = set(v[0] for v in cursor.fetchall())
            self.assertFalse(any(v in keys for v in ('user', 'password')))

    def testNoSecretsInConfig(self):
        """Make sure configure_from_database isn't retrieveing accidentaly stored secrets."""
        utils._initdb(config=self._config, confirm=True)
        with utils.connect(self._config).cursor() as cursor:
            for key in ('password', 'user'):
                cursor.execute(
                    'insert into registration_metadata values (%s,%s);',
                    (key, self._config[key]))
        config_from_database = utils.configure_from_database(
            dbname=self._config['dbname'], dbtype=self._config['dbtype'])
        self.assertFalse(
            any(v in config_from_database for v in ('user', 'password')))


if __name__ == '__main__':
    unittest.main()
