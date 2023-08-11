# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest
import sqlite3
from rdkit import Chem
from rdkit.Chem import rdDistGeom
import random
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
    cfg['host'] = 'localhost'
    cfg['dbtype'] = 'postgresql'
    try:
        cn = utils._connect(config=cfg)
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
        utils.initdb(config=self._config, confirm=True)
        for smi in smis:
            utils.register(smiles=smi, config=self._config)
        mols = [Chem.MolFromSmiles(x) for x in ('Cc1[nH]ncc1', 'Cc1n[nH]cc1')]
        for mol in mols:
            utils.register(mol=mol, config=self._config)

    def testRegister(self):
        utils.initdb(config=self._config, confirm=True)
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
        utils.initdb(config=self._config, confirm=True)
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

    def testBulkRegisterAllowDupes(self):
        utils.initdb(config=self._config, confirm=True)
        mols = [
            Chem.MolFromSmiles(x)
            for x in ('CCC', 'CCCO', 'C1', 'c1cc1', 'CCC', 'C1CC1')
        ]

        res = utils.bulk_register(mols=mols,
                                  failOnDuplicate=False,
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
        tpl = res[0]
        self.assertEqual(len(tpl), 3)
        self.assertEqual(tpl[0], 1)
        mb = Chem.MolFromMolBlock(tpl[1])
        self.assertEqual(
            utils.query(smiles=Chem.MolToSmiles(mb), config=self._config), [1])
        self.assertEqual(tpl[2], 'mol')
        res = utils.retrieve(ids=[100], config=self._config)
        self.assertEqual(len(res), 0)

    def testStandardizationOptions(self):
        lconfig = self._config.copy()
        lconfig['standardization'] = 'charge'
        utils.initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='CCCO', config=lconfig), 1)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='CCC[O-]', config=lconfig))

        lconfig['standardization'] = 'tautomer'
        utils.initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='Cc1[nH]ncc1', config=lconfig),
                         1)
        self.assertRaises(
            self.integrityError,
            lambda: utils.register(smiles='Cc1n[nH]cc1', config=lconfig))

    def testStandardizationFunctions(self):
        lconfig = self._config.copy()

        def evenAtomCount(mol):
            if mol.GetNumAtoms() % 2:
                return None
            return mol

        # Silly example which only accepts even numbers of atoms
        lconfig['standardization'] = evenAtomCount

        utils.initdb(config=lconfig, confirm=True)
        self.assertEqual(utils.register(smiles='CCCO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='CCC', config=lconfig),
                         RegistrationFailureReasons.FILTERED)

        utils.initdb(config=lconfig, confirm=True)
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
        utils.initdb(config=lconfig, confirm=True)
        lconfig['standardization'] = [evenAtomCount, hasAnOxygen]
        self.assertEqual(utils.register(smiles='CC', config=lconfig),
                         RegistrationFailureReasons.FILTERED)
        self.assertEqual(utils.register(smiles='CO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='CCO', config=lconfig),
                         RegistrationFailureReasons.FILTERED)

        # ----------------------
        # mixing standardization options and functions
        utils.initdb(config=lconfig, confirm=True)
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

        utils.initdb(config=lconfig, confirm=True)
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


class TestLWRegTautomerv2(unittest.TestCase):
    integrityError = sqlite3.IntegrityError

    def setUp(self):
        cn = sqlite3.connect(':memory:')
        self._config = utils.defaultConfig()
        self._config['connection'] = cn
        self._config['useTautomerHashv2'] = 1

    def testRegister(self):
        utils.initdb(config=self._config, confirm=True)
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


class TestStandardizationLabels(unittest.TestCase):
    def testStandards(self):
        cfg = utils.defaultConfig()
        for k in utils.standardizationOptions:
            cfg['standardization'] = k
            lbl = utils._get_standardization_label(cfg)
            self.assertEqual(lbl, k)

    def testCombined(self):
        cfg = utils.defaultConfig()
        for k in utils.standardizationOptions:
            cl = ['foo', k]
            cfg['standardization'] = cl
            lbl = utils._get_standardization_label(cfg)
            self.assertEqual(lbl, '|'.join(cl))

    def testOthers(self):
        def func1(x):
            pass

        cfg = utils.defaultConfig()
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
        cfg = utils.defaultConfig()
        cfg['dbname'] = 'foo.sqlt'
        utils.initdb(config=cfg, confirm=True)
        self.assertEqual(utils.register(smiles='CCO', config=cfg), 1)
        self.assertEqual(utils.register(smiles='CCOC', config=cfg), 2)
        oac = standardization_lib.OverlappingAtomsCheck()
        cfg['standardization'] = ['fragment', oac]
        self.assertEqual(utils.register(smiles='CCN', config=cfg), 3)
        self.assertEqual(utils.register(smiles='CCNC', config=cfg), 4)
        cn = utils._connect(cfg)
        curs = cn.cursor()
        curs.execute(
            "select count(*) from molblocks where standardization is not null")
        self.assertEqual(curs.fetchone()[0], 2)
        curs.execute(
            "select count(*) from molblocks where standardization is null")
        self.assertEqual(curs.fetchone()[0], 2)


class TestConformerHashes(unittest.TestCase):
    def setUp(self):
        self._config = utils.defaultConfig()
        self._config['hashConformer'] = True
        self._mol1 = Chem.AddHs(Chem.MolFromSmiles('OC(=O)CCCC'))
        rdDistGeom.EmbedMolecule(self._mol1, randomSeed=0xf00d)
        self._mol2 = Chem.Mol(self._mol1)
        rdDistGeom.EmbedMolecule(self._mol2, randomSeed=0xf00d + 1)

    def testConformerDupes(self):
        utils.initdb(config=self._config, confirm=True)
        self.assertEqual(utils.register(mol=self._mol1, config=self._config),
                         1)
        self.assertEqual(utils.register(mol=self._mol2, config=self._config),
                         2)
        aorder = list(range(self._mol1.GetNumAtoms()))
        random.shuffle(aorder)
        nmol = Chem.RenumberAtoms(self._mol1, aorder)
        self.assertEqual(
            utils.register(mol=nmol,
                           config=self._config,
                           fail_on_duplicate=False), 1)


if __name__ == '__main__':
    unittest.main()