# LWReg: a lightweight chemical registration system

This provides a basic registration system which can be used either as a python
library or via a command-line interface.

Basic operations:
- `initdb`: resets the database. **Note** that this destroys any information which is already in the database, so be careful with it.
- `register`: standardize the input molecule, calculates a hash for it, and adds them molecule to the database if it's not already there. returns the molregno (registry ID) of the newly registered molecule.
- `query`: takes a molecule as input and checks whether or not a matching molecule is registered. returns molregnos (registry IDs) of the matching molecule(s), if any
- `retrieve`: takes one or more IDs and returns the registered structures for them

## Dependencies:
- rdkit
- click
- psycopg2 (only if you want to use a postgresql database)

## Install:
After checking out this repo, run this command in this directory:
```
pip install --editable .
```

## Very basic usage demo

### Command line
```
% lwreg initdb --confirm=yes
% lwreg register --smiles CCOCC
1
% lwreg register --smiles CCOCCC
2
% lwreg register --smiles CCNCCC
3
% lwreg register --smiles CCOCCC
Traceback (most recent call last):
  File "/home/glandrum/miniconda3/envs/py310_rdkit/bin/lwreg", line 33, in <module>
    sys.exit(load_entry_point('lwreg', 'console_scripts', 'lwreg')())
  File "/home/glandrum/miniconda3/envs/py310_rdkit/lib/python3.10/site-packages/click/core.py", line 1128, in __call__
    return self.main(*args, **kwargs)
  File "/home/glandrum/miniconda3/envs/py310_rdkit/lib/python3.10/site-packages/click/core.py", line 1053, in main
    rv = self.invoke(ctx)
  File "/home/glandrum/miniconda3/envs/py310_rdkit/lib/python3.10/site-packages/click/core.py", line 1659, in invoke
    return _process_result(sub_ctx.command.invoke(sub_ctx))
  File "/home/glandrum/miniconda3/envs/py310_rdkit/lib/python3.10/site-packages/click/core.py", line 1395, in invoke
    return ctx.invoke(self.callback, **ctx.params)
  File "/home/glandrum/miniconda3/envs/py310_rdkit/lib/python3.10/site-packages/click/core.py", line 754, in invoke
    return __callback(*args, **kwargs)
  File "/home/glandrum/Code/lightweight-registration/lwreg/lwreg.py", line 78, in register
    return utils.register(**kwargs)
  File "/home/glandrum/Code/lightweight-registration/lwreg/utils.py", line 249, in register
    mrn = _register_mol(tpl, escape, cn, curs, config, fail_on_duplicate)
  File "/home/glandrum/Code/lightweight-registration/lwreg/utils.py", line 184, in _register_mol
    curs.execute(
sqlite3.IntegrityError: UNIQUE constraint failed: hashes.fullhash
% lwreg query --smiles CCOCCC
2
% lwreg retrieve --id 2
(2, '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 6 5 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 5 C 5.196152 -0.000000 0.000000 0\nM  V30 6 C 6.495191 0.750000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 4 5\nM  V30 5 1 5 6\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n', 'mol')
```

### Python
```
In [1]: import lwreg

In [3]: lwreg.initdb(confirm=True)
Out[3]: True

In [4]: lwreg.register(smiles='CCO')
Out[4]: 1

In [5]: lwreg.register(smiles='CCOC')
Out[5]: 2

In [6]: from rdkit import Chem

In [7]: m = Chem.MolFromSmiles('CCOCC')

In [8]: lwreg.register(mol=m)
Out[8]: 3

In [9]: lwreg.register(mol=m)
---------------------------------------------------------------------------
IntegrityError                            Traceback (most recent call last)
Input In [9], in <cell line: 1>()
----> 1 lwreg.register(mol=m)

File ~/Code/lightweight-registration/lwreg/utils.py:249, in register(config, mol, molfile, molblock, smiles, escape, fail_on_duplicate, no_verbose)
    247 cn = _connect(config)
    248 curs = cn.cursor()
--> 249 mrn = _register_mol(tpl, escape, cn, curs, config, fail_on_duplicate)
    250 if not no_verbose:
    251     print(mrn)

File ~/Code/lightweight-registration/lwreg/utils.py:184, in _register_mol(tpl, escape, cn, curs, config, failOnDuplicate)
    181     mhash, layers = hash_mol(sMol, escape=escape, config=config)
    183     # will fail if the fullhash is already there
--> 184     curs.execute(
    185         _replace_placeholders(
    186             'insert into hashes values (?,?,?,?,?,?,?,?,?)'), (
    187                 mrn,
    188                 mhash,
    189                 layers[RegistrationHash.HashLayer.FORMULA],
    190                 layers[RegistrationHash.HashLayer.CANONICAL_SMILES],
    191                 layers[RegistrationHash.HashLayer.NO_STEREO_SMILES],
    192                 layers[RegistrationHash.HashLayer.TAUTOMER_HASH],
    193                 layers[RegistrationHash.HashLayer.NO_STEREO_TAUTOMER_HASH],
    194                 layers[RegistrationHash.HashLayer.ESCAPE],
    195                 layers[RegistrationHash.HashLayer.SGROUP_DATA],
    196             ))
    198     cn.commit()
    199 except _violations:

IntegrityError: UNIQUE constraint failed: hashes.fullhash

In [10]: lwreg.query(smiles='CCOC')
Out[10]: [2]

In [11]: lwreg.query(smiles='CCOCC')
Out[11]: [3]

In [12]: lwreg.query(smiles='CCOCO')
Out[12]: []

In [13]: lwreg.retrieve(id=2)
Out[13]: 
((2,
  '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n',
  'mol'),)

In [14]: lwreg.retrieve(ids=[2,3])
Out[14]: 
((2,
  '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n',
  'mol'),
 (3,
  '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 5 4 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 5 C 5.196152 -0.000000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 4 5\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n',
  'mol'))


```
