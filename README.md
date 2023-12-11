# LWReg: a lightweight chemical registration system

This provides a basic registration system which can be used either as a python
library or via a command-line interface.

Basic operations:
- `initdb`: resets the database. **Note** that this destroys any information which is already in the database, so be careful with it.
- `register`: standardize the input molecule, calculates a hash for it, and adds them molecule to the database if it's not already there. returns the molregno (registry ID) of the newly registered molecule.
- `query`: takes a molecule as input and checks whether or not a matching molecule is registered. returns molregnos (registry IDs) of the matching molecule(s), if any
- `retrieve`: takes one or more IDs and returns the registered structures for them

## Installation for non-experts:

Please look at the INSTALL.md file.

## Dependencies:
- rdkit v2023.03.1 or later
- click
- psycopg2 (only if you want to use a postgresql database)

## Installation for experts:
After installing the dependencies (above) and checking out this repo, run this command in this directory:
```
pip install --editable .
```

## Run in Docker
```shell
docker build -t lwreg .

# Run Jupyter notebook on the docker container
docker run -i -t -p 8888:8888 rdkit-lwreg /bin/bash -c "\
    apt update && apt install libtiff5 -y && \
    pip install notebook && \
    jupyter notebook \
    --notebook-dir=/lw-reg --ip='*' --port=8888 \
    --no-browser --allow-root"
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
ERROR:root:Compound already registered
% lwreg query --smiles CCOCCC
2
% lwreg retrieve --id 2
(2, '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 6 5 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 O 2.598076 -0.000000 0.000000 0\nM  V30 4 C 3.897114 0.750000 0.000000 0\nM  V30 5 C 5.196152 -0.000000 0.000000 0\nM  V30 6 C 6.495191 0.750000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 4 5\nM  V30 5 1 5 6\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n', 'mol')
```

### Python
```
In [1]: import lwreg

In [3]: lwreg.initdb()
This will destroy any existing information in the registration database.
  are you sure? [yes/no]: yes
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

## Custom standardization/filtering functions

When using the Python API you have extensive control over the standardization and validation operations which are performed on the molecule.

Start with a couple of examples showing what the 'fragment' and 'charge' built-in standardizers do:

```
>>> config['standardization'] = 'fragment'
>>> Chem.MolToSmiles(lwreg.utils.standardize_mol(Chem.MolFromSmiles('CC[O-].[Na+]'),config=config))
'CC[O-]'
>>> config['standardization'] = 'charge'
>>> Chem.MolToSmiles(lwreg.utils.standardize_mol(Chem.MolFromSmiles('CC[O-].[Na+]'),config=config))
'CCO'
```

Now define a custom filter which rejects (by returning None) molecules which have a net charge and then use that:
```
>>> def reject_charged_molecules(mol):
...     if Chem.GetFormalCharge(mol):
...         return None
...     return mol
...
>>> config['standardization'] = reject_charged_molecules
>>> Chem.MolToSmiles(lwreg.utils.standardize_mol(Chem.MolFromSmiles('CC[O-].[Na+]'),config=config))
'CC[O-].[Na+]'
```

Here's an example which fails:

```
>>> lwreg.utils.standardize_mol(Chem.MolFromSmiles('CC[O-]'),config=config) is None
True
```

We can chain standardization/filtering operations together by providing a list. The individual operations are run in order. Here's an example where we attempt to neutralise the molecule by finding the charge parent and then apply our reject_charged_molecules filter:
```
>>> config['standardization'] = ['charge',reject_charged_molecules]
>>> lwreg.utils.standardize_mol(Chem.MolFromSmiles('CC[O-]'),config=config) is None
False
>>> lwreg.utils.standardize_mol(Chem.MolFromSmiles('CC[N+](C)(C)C'),config=config) is None
True
```
That last one failed because the quarternary nitrogen can't be neutralized.

There are a collection of other standardizers/filters available in the module lwreg.standardization_lib


## Registering conformers

When the configuration option `registerConformers` is set to True, lwreg expects that the compounds to be registered will have an associated conformer. The conformers are tracked in a different table than the molecule topologies and expectation is that every molecule registered will have a conformer (it's an error if they don't). It is possible to register multiple conformers for a single molecular structure (topology).

Note that once a database is created in `registerConformers` mode, it probably should always be used in that mode. 

### Differences when in `registerConformers` mode

- `register()` and `bulk_register()` require molecules to have associated conformers. Both return `(molregno, conf_id)` tuples instead of just `molregno`s
- `query()`: if called with the `ids` argument, this will return all of the conformers for the supplied molregnos as `(molregno, conf_id)` tuples. If called with a molecule, the conformer of the molecule will be hashed and looked up in the `conformers`` table, returns a list of `(molregno, conf_id)` tuples.
- `retrieve()`: if called with `(molregno, conf_id)` tuple(s), this will return `(molregno, conf_id, molblock)` tuples where the `molblock`s contain the coordinates of the registered conformers.

### Hashing conformers for registration

Just as molecular hashes are used to recognize when two molecules are the same, lwreg uses a hashing scheme to detect when two conformers are the same. The algorithm for this is simple:
The atomic positions are converted into strings (rounding the floating point values to a fixed, but configurable, number of digits), sorting the positions, and then combining them into a single string, which is the final hash.

If registering a multi-conformer molecule, it is most efficient to call `register_multiple_conformers()`. That only does the work of standardizing the molecule and calculating the molecule hash once.

# Data layout

## The base tables
Here's the SQL to create the base lwreg tables in sqlite:
```
create table registration_metadata (key text, value text);
create table hashes (molregno integer primary key, fullhash text unique, 
            formula text, canonical_smiles text, no_stereo_smiles text, 
            tautomer_hash text, no_stereo_tautomer_hash text, "escape" text, sgroup_data text, rdkitVersion text);
create table orig_data (molregno integer primary key, data text, datatype text);
create table molblocks (molregno integer primary key, molblock text, standardization text);
```


## The conformers table
Here's the SQL to create the conformers table in sqlite when `registerConformers` is set:
```
create table conformers (conf_id integer primary key, molregno integer not null, 
                   conformer_hash text not null unique, molblock text);
```