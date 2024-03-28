# LWReg: a lightweight chemical registration system

This provides a basic registration system which can be used either as a python
library or via a command-line interface.

Basic operations:
- `initdb`: resets the database. **Note** that this destroys any information which is already in the database, so be careful with it.
- `register`: standardize the input molecule, calculates a hash for it, and adds them molecule to the database if it's not already there. returns the molregno (registry ID) of the newly registered molecule.
- `query`: takes a molecule as input and checks whether or not a matching molecule is registered. returns molregnos (registry IDs) of the matching molecule(s), if any
- `retrieve`: takes one or more IDs and returns the registered structures for them

## Quick-start installation for non-experts:

Assuming that you have conda (or mamba or something equivalent) installed you can install lwreg directly from this github repo by first creating a conda environment with all the dependencies installed:
```
% conda env create --name py311_lwreg --file=https://raw.githubusercontent.com/rinikerlab/lightweight-registration/main/environment.yml
```
If you have mamba installed, you can run this instead (it will run faster):
```
% mamba env create --name py311_lwreg --file=https://raw.githubusercontent.com/rinikerlab/lightweight-registration/main/environment.yml
```

You can then activate the new environment and install lwreg:
```
% conda activate py311_lwreg
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

## Installation for non-experts:

Please look at the INSTALL.md file.

## Installation for experts:

### Dependencies:
- rdkit v2023.03.1 or later
- click
- psycopg2 (only if you want to use a postgresql database)

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
Using the default configuration options:
```
% lwreg initdb --confirm=yes                          # initialization with default config
```
Instead using the interactive configuration:
```
% lwreg interactive-config                            # interactive configuration, add host, user, password to config.json afterwards
% lwreg --config config.json initdb --confirm=yes     # initialization with config.json
```
```
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
>>> import lwreg

>>> from lwreg import utils, helpers
```
Using the default configuration options:
```
# default config 
>>> lwreg.set_default_config(utils.defaultConfig())   # you generally will want to provide more information about the database

>>> lwreg.initdb()
This will destroy any existing information in the registration database.
  are you sure? [yes/no]: yes
True
```
Instead using the interactive configuration:
```
>>> config = helpers.interactive_config()   # interactive configuration

>>> lwreg.initdb(config=config)
This will destroy any existing information in the registration database.
  are you sure? [yes/no]: yes
True
```
```
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

## Custom standardization/filtering functions

When using the Python API you have extensive control over the standardization and validation operations which are performed on the molecule.

Start with a couple of examples showing what the 'fragment' and 'charge' built-in standardizers do:

```
>>> config['standardization'] = 'fragment'
>>> Chem.MolToSmiles(lwreg.standardize_mol(Chem.MolFromSmiles('CC[O-].[Na+]'),config=config))
'CC[O-]'
>>> config['standardization'] = 'charge'
>>> Chem.MolToSmiles(lwreg.standardize_mol(Chem.MolFromSmiles('CC[O-].[Na+]'),config=config))
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
>>> Chem.MolToSmiles(lwreg.standardize_mol(Chem.MolFromSmiles('CC[O-].[Na+]'),config=config))
'CC[O-].[Na+]'
```

Here's an example which fails:

```
>>> lwreg.standardize_mol(Chem.MolFromSmiles('CC[O-]'),config=config) is None
True
```

We can chain standardization/filtering operations together by providing a list. The individual operations are run in order. Here's an example where we attempt to neutralise the molecule by finding the charge parent and then apply our reject_charged_molecules filter:
```
>>> config['standardization'] = ['charge',reject_charged_molecules]
>>> lwreg.standardize_mol(Chem.MolFromSmiles('CC[O-]'),config=config) is None
False
>>> lwreg.standardize_mol(Chem.MolFromSmiles('CC[N+](C)(C)C'),config=config) is None
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
- `retrieve()`: if called with `(molregno, conf_id)` tuple(s), this will return a dictionary of `(molblock, 'mol')` tuples with `(molregno, conf_id)` tuples as keys where the `molblock`s contain the coordinates of the registered conformers.

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
create table orig_data (molregno integer primary key, data text, datatype text, timestamp DATETIME DEFAULT CURRENT_TIMESTAMP);
create table molblocks (molregno integer primary key, molblock text, standardization text);
```


## The conformers table
Here's the SQL to create the conformers table in sqlite when `registerConformers` is set:
```
create table conformers (conf_id integer primary key, molregno integer not null, 
                   conformer_hash text not null unique, molblock text);
```