# LWReg: a lightweight chemical registration system

[![Lwreg](https://img.shields.io/badge/DOI-10.1021/acs.jcim.3c00800-blue)](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01133)
[![ReadTheDocs](https://readthedocs.org/projects/docs/badge/?version=latest)](https://lightweight-registration.readthedocs.io/en/)
[![RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
![TOC Image](https://github.com/rinikerlab/lightweight-registration/blob/feat/docs/Lwreg.png)

This provides a basic registration system which can be used either as a python
library or via a command-line interface.

Basic operations:
- `initdb`: resets the database. **Note** that this destroys any information which is already in the database, so be careful with it.
- `register`: standardize the input molecule, calculates a hash for it, and adds them molecule to the database if it's not already there. returns the molregno (registry ID) of the newly registered molecule.
- `query`: takes a molecule as input and checks whether or not a matching molecule is registered. returns molregnos (registry IDs) of the matching molecule(s), if any
- `retrieve`: takes one or more IDs and returns the registered structures for them

## Publications

[1] J. Chem. Inf. Model. 2024, 64, 16, 6247–6252 : https://pubs.acs.org/doi/10.1021/acs.jcim.4c01133

## Documentation

Our in detail documentation: https://lightweight-registration.readthedocs.io/en/

## Quick-start

Assuming that you have conda (or mamba or something equivalent) installed you can install lwreg directly from this github repo by first creating a conda environment with all the dependencies installed:
```
% conda env create --name py313_lwreg --file=https://raw.githubusercontent.com/rinikerlab/lightweight-registration/main/environment.yml
```
If you have mamba installed, you can run this instead (it will run faster):
```
% mamba env create --name py313_lwreg --file=https://raw.githubusercontent.com/rinikerlab/lightweight-registration/main/environment.yml
```

You can then activate the new environment and install lwreg:
```
% conda activate py313_lwreg
% python -m pip install git+https://github.com/rinikerlab/lightweight-registration
```

You can then verify that the install worked by doing:
```
% lwreg --help
```

If you want to use PostgreSQL as the database for lwreg, then you will also need to install the python connector for PostgreSQL:
```
% conda install -c conda-forge psycopg2
```

For further information, consult the INSTALL.md file.

## Very basic usage demo

### Command line
```
% lwreg initdb --confirm=yes
% lwreg register --smiles CCOCC
% lwreg query --smiles CCOCCC
% lwreg retrieve --id 2
```

### Python
```
>>> import lwreg

>>> from lwreg import utils

>>> lwreg.set_default_config(utils.defaultConfig())   # you generally will want to provide more information about the database

>>> lwreg.initdb()
This will destroy any existing information in the registration database.
  are you sure? [yes/no]: yes
True

>>> lwreg.register(smiles='CCO')
1

>>> lwreg.register(smiles='CCOC')
2

>>> from rdkit import Chem

>>> m = Chem.MolFromSmiles('CCOCC')

>>> lwreg.register(mol=m)
3

>>> lwreg.register(mol=m)
---------------------------------------------------------------------------
IntegrityError                            Traceback (most recent call last)
Input In [10], in <cell line: 1>()
----> 1 lwreg.register(mol=m)
  
  ... DETAILS REMOVED ...

IntegrityError: UNIQUE constraint failed: hashes.fullhash

>>> lwreg.query(smiles='CCOC')
[2]

>>> lwreg.query(smiles='CCOCC')
[3]

>>> lwreg.query(smiles='CCOCO')
[]

>>> lwreg.retrieve(id=2)
{2: ('\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n',
  'mol')}

>>> lwreg.retrieve(ids=[2,3])
{2: ('\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n',
  'mol'),
 3: ('\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 5 4 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 5 C 5.196152 -0.000000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 4 5\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n',
  'mol')}


```