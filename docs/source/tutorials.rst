Tutorials
=========

Here, we provide a small set of tutorials to help you get started with lwreg. 
Please keep in mind that lwreg allows for problem dependent customization:
The tutorials are meant as a starting point, but might not be directly applicable to your specific use case.

Setting up a Database and First Registration
---------------------------------------------
For a single user, SQLite is the easiest way to get started with lwreg.
In this first example, we will show how to set up an lwreg instance for registering molecules and querying them.
We will also see how to retrieve the registerd structures::
     
  import lwreg
  from lwreg import utils

  # we need to start with a config, here let's use the default provided by lwreg
  config = utils.defaultConfig()

  # the default config is the following:
  # config = {'dbname': './testdb.sqlt',
  # 'dbtype': 'sqlite3',
  # 'standardization': 'fragment',
  # 'removeHs': 1,
  # 'useTautomerHashv2': 0,
  # 'registerConformers': 0,
  # 'numConformerDigits': 3,
  # 'lwregSchema': ''}

  # now we can initialize the database
  lwreg.initdb(config)

The user will need to confirm the initialization of the database. 
After that, the database is ready to be used and we can proceed to register molecules. ::

  # after confirming the intialization, we can register some molecules
  from rdkit import Chem
  from rdkit.Chem import Draw
  from rdkit.Chem.Draw import IPythonConsole
  
  smis = ["Cc1[nH]ncc1[C@H](F)Cl","Cc1[nH]ncc1[C@@H](F)Cl","Cc1[nH]ncc1[CH](F)Cl","Cc1n[nH]cc1[C@@H](F)Cl"]
  ms = [Chem.MolFromSmiles(smi) for smi in smis]
  IPythonConsole.drawOptions.useBWAtomPalette()
  IPythonConsole.drawOptions.legendFontSize=24
  Draw.MolsToGridImage(ms,legends=[str(i+1) for i in range(len(ms))],molsPerRow=2,subImgSize=(300,300))

You will see the molecules drawn on a grid displayed in the notebook.
Let's also assume that you set up your database and are now coming back to register some molecules.
In that case, you should always start by retrieving the config from the database. ::

  config = utils.configure_from_database(dbname='./testdb.sqlt',dbtype='sqlite3')

  # registering one molecule
  lwreg.register(config,ms[0])
  # registering multiple molecules
  lwreg.bulk_register(config,ms[1::])

  # there should now be 4 molecules in the database
  utils.get_all_registry_numbers(config)

  # we can now query the database
  lwreg.query(smiles='Cc1[nH]ncc1[C@H](F)Cl')

  # query using just the no_stereo_smiles layer (this pays attention to tautomers but ignores stereochemistry):
  ids = lwreg.query(smiles='Cc1[nH]ncc1[C@H](F)Cl',
            layers=utils.HashLayer.NO_STEREO_SMILES)

  # retrieving the structures for the query matches
  res = lwreg.retrieve(ids=ids)

  # retrieve returns a dictionary with the molregnos as keys. The values are two-tuples with the molecule structure and its configure_from_database
  res[1] 

  # build an rdkit mol object from the retrieved information
  mol = Chem.MolFromMolBlock(res[1][0],removeHs=False)


Tutorial Notebooks
------------------

* `01_registration_basics.ipynb <https://github.com/rinikerlab/lightweight-registration/blob/main/tutorial/01_registration_basics.ipynb>`_
   Setting up a database and the basics of registering and querying with molecules.
* `02_querying_existing_database.ipynb <https://github.com/rinikerlab/lightweight-registration/blob/main/tutorial/02_querying_existing_database.ipynb>`_
   Connecting to an existing database and querying for and retrieving registered molecules._
* `03_registering_conformers.ipynb <https://github.com/rinikerlab/lightweight-registration/blob/main/tutorial/03_registering_conformers.ipynb>`_
   Working with lwreg in conformer mode to store conformers with the molecules.
* `04_registering_molecules_and_data.ipynb <https://github.com/rinikerlab/lightweight-registration/blob/main/tutorial/04_registering_molecules_and_data.ipynb>`_
   Storing additional data together with the registered molecules.

Advanced Tutorials
^^^^^^^^^^^^^^^^^^

* `100_computational_experiment_1.ipynb <https://github.com/rinikerlab/lightweight-registration/blob/main/tutorial/100_computational_experiment_1.ipynb>`_
   Storing and working with results from a computational experiment that produces files.
* `101_computational_experiment_2.ipynb <https://github.com/rinikerlab/lightweight-registration/blob/main/tutorial/101_computational_experiment_2.ipynb>`_
   Storing and working with results from a computational experiment that uses a multi-conformer database.