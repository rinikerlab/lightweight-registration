# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest
import sqlite3
from rdkit import Chem
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
    cfg['dbname'] = 'dbname=lwreg_tests host=localhost'
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
        self.assertEqual(utils.register(smiles='CCCOC', config=self._config),
                         3)
        self.assertEqual(
            utils.register(mol=Chem.MolFromSmiles('C1CCC1'),
                           config=self._config), 4)
        self.assertEqual(utils.register(mol=None, config=self._config),
                         RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(utils.register(smiles='CCCC1', config=self._config),
                         RegistrationFailureReasons.PARSE_FAILURE)
        self.assertEqual(utils.register(smiles='c1nccc1', config=self._config),
                         RegistrationFailureReasons.PARSE_FAILURE)

    def testBulkRegister(self):
        utils.initdb(config=self._config, confirm=True)
        mols = [
            Chem.MolFromSmiles(x)
            for x in ('CCC', 'CCCO', 'C1', 'c1cc1', 'CCC', 'C1CC1')
        ]
        self.assertEqual(utils.bulk_register(mols=mols, config=self._config),
                         (1, 2, RegistrationFailureReasons.PARSE_FAILURE,
                          RegistrationFailureReasons.PARSE_FAILURE,
                          RegistrationFailureReasons.DUPLICATE, 3))

    def testBulkRegisterAllowDupes(self):
        utils.initdb(config=self._config, confirm=True)
        mols = [
            Chem.MolFromSmiles(x)
            for x in ('CCC', 'CCCO', 'C1', 'c1cc1', 'CCC', 'C1CC1')
        ]
        self.assertEqual(
            utils.bulk_register(mols=mols,
                                failOnDuplicate=False,
                                config=self._config),
            (1, 2, RegistrationFailureReasons.PARSE_FAILURE,
             RegistrationFailureReasons.PARSE_FAILURE, 1, 3))

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
        self.assertEqual(tpl[0], 1)
        mb = Chem.MolFromMolBlock(tpl[1])
        self.assertEqual(
            utils.query(smiles=Chem.MolToSmiles(mb), config=self._config), [1])
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
        self.assertIsNone(utils.register(smiles='CCC', config=lconfig))

        utils.initdb(config=lconfig, confirm=True)
        checker = standardization_lib.OverlappingAtomsCheck()
        lconfig['standardization'] = checker
        self.assertEqual(
            utils.register(smiles='CC |(0,0,;0,1,)|', config=lconfig), 1)
        self.assertEqual(
            utils.register(smiles='CO |(0,0,;0,0,)|', config=lconfig), None)

        def hasAnOxygen(mol):
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 8:
                    return mol
            return None

        # ----------------------
        # chaining standardization functions
        utils.initdb(config=lconfig, confirm=True)
        lconfig['standardization'] = [evenAtomCount, hasAnOxygen]
        self.assertEqual(utils.register(smiles='CC', config=lconfig), None)
        self.assertEqual(utils.register(smiles='CO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='CCO', config=lconfig), None)

        # ----------------------
        # mixing standardization options and functions
        utils.initdb(config=lconfig, confirm=True)
        lconfig['standardization'] = ('fragment', evenAtomCount, hasAnOxygen)
        self.assertEqual(utils.register(smiles='CC', config=lconfig), None)
        self.assertEqual(utils.register(smiles='CO', config=lconfig), 1)
        self.assertEqual(utils.register(smiles='C[O-].[Na]', config=lconfig),
                         2)
        self.assertEqual(utils.register(smiles='CCO', config=lconfig), None)

    def testStandardizationLibCheckers(self):
        lconfig = self._config.copy()

        utils.initdb(config=lconfig, confirm=True)
        checker = standardization_lib.OverlappingAtomsCheck()
        lconfig['standardization'] = checker
        self.assertEqual(
            utils.register(smiles='CC |(0,0,;0,1,)|', config=lconfig), 1)
        self.assertEqual(
            utils.register(smiles='CO |(0,0,;0,0,)|', config=lconfig), None)

        checker = standardization_lib.PolymerCheck()
        lconfig['standardization'] = checker
        self.assertEqual(
            utils.register(smiles='*-CO-* |$star_e;;;star_e$,Sg:n:1,2::ht|',
                           config=lconfig), None)
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
        self._config['dbname'] = 'dbname=lwreg_tests host=localhost'
        self._config['dbtype'] = 'postgresql'


if __name__ == '__main__':
    unittest.main()