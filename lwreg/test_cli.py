# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import unittest
from click.testing import CliRunner
from . import lwreg
import tempfile
import os


# just some basic testing to make sure that the CLI works.
# The real testing of the underlying library is in test_lwreg.py
class TestLWRegCLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        jsonf = tempfile.NamedTemporaryFile(suffix='.json',
                                            mode='w+t',
                                            delete=False)
        sqltf = tempfile.NamedTemporaryFile(suffix='.sqlt', delete=False)
        cls.sqltf = sqltf.name.replace('\\', '/')

        jsonf.write('''{"dbname": "%s"}''' % (cls.sqltf, ))
        jsonf.flush()
        cls.configFile = jsonf.name
        print(cls.configFile)
        jsonf.close()

    @classmethod
    def __tearDownClass(cls):
        try:
            os.unlink(cls.configFile)
        except:
            pass
        try:
            os.unlink(cls.sqltf)
        except:
            pass

    def test1(self):
        runner = CliRunner()
        result = runner.invoke(
            lwreg.cli,
            [f'--config={self.configFile}', 'initdb', '--confirm=yes'])
        self.assertEqual(result.exit_code, 0)
        result = runner.invoke(
            lwreg.cli,
            [f'--config={self.configFile}', 'register', '--smiles=CCC'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), '1')

        # expected failure
        result = runner.invoke(
            lwreg.cli,
            [f'--config={self.configFile}', 'register', '--smiles=CCC'])
        self.assertEqual(result.exit_code, 1)
        self.assertEqual(result.output.strip(), '')

        result = runner.invoke(lwreg.cli, [
            f'--config={self.configFile}', 'register', '--smiles=CCC',
            '--fail-on-duplicate'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), '1')

        result = runner.invoke(
            lwreg.cli,
            [f'--config={self.configFile}', 'register', '--smiles=CCCO'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), '2')

        result = runner.invoke(lwreg.cli, [
            f'--config={self.configFile}', 'register', '--smiles=CCC',
            '--escape=really'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), '3')

        result = runner.invoke(lwreg.cli, [
            f'--config={self.configFile}',
            'query',
            '--smiles=COC',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), 'not found')

        result = runner.invoke(lwreg.cli, [
            f'--config={self.configFile}',
            'query',
            '--smiles=CCC',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), '1')

        result = runner.invoke(lwreg.cli, [
            f'--config={self.configFile}', 'query', '--smiles=CCC',
            '--layers=FORMULA'
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), '1 3')

        result = runner.invoke(lwreg.cli, [
            f'--config={self.configFile}',
            'retrieve',
            '--id=1',
        ])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("(1, ", result.output)
        self.assertIn("M  END\\n',", result.output)
        self.assertIn("'mol'", result.output)


if __name__ == '__main__':
    unittest.main()