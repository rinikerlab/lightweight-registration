Introduction
=============

.. _GetStarted:


What is a registration system?
------------------------------
A chemical registration system is a software system (e.g. database) that holds chemical structures and assigns them identifiers.
In order to do so, a chemical registration system needs to be able to check a presented compound against already registered ones and decide upon its uniqueness. 
The identifiers can be used as keys when storing additional data about the chemical structures in additional tables.

lwreg's Basic Functionality
---------------------------
There are four basic operations provides by the lwreg package that can either by called via the Python API or the command line interface:

.. list-table:: Basic operations
   :widths: 10 30
   :header-rows: 1

   * - Operation
     - Description
   * - initdb
     - Initialize a new database. **Warning**: This will delete all existing data in the database.
   * - register
     - Register a new compound, returns a unique identifier (molregno) upon success.
   * - query
     - Query the database if a compound is already registered.
       If so, returns the respective molregno(s) matching the search.
   * - retrieve
     - Get the resgisted structure(s) for one or multiple given identifier(s). Ids in the retrieve list that do not exist do not show up in the results.

See the :ref:`Basic Functionality Reference` for more details and :ref:`Additional Functionality Reference` for advanced features.

Installation Guide
------------------
The following outline assumes that the user has conda installed. 

.. code-block:: bash

    conda env create --name py313_lwreg --file=https://raw.githubusercontent.com/rinikerlab/lightweight-registration/main/environment.yml
    conda activate py313_lwreg
    python -m pip install git+https://github.com/rinikerlab/lightweight-registration

To verify the installation by running

.. code-block:: bash

    lwreg --help

When starting to use lwreg, you also must choose your prefered database management system. 
The two systems currently supported by lwreg are SQLite and PostgreSQL.

SQLite is a serverless database system that stores the database in a single file and which is supported as part of the Python standard library.

If you choose to use PostgreSQL, you need to install the psycopg package (note that this is done automatically if you use the environment setup described above):

.. code-block:: bash

    conda install -c conda-forge psycopg

You will also need to have the PostgreSQL server installed and running.
Please refer to the `PostgreSQL documentation <https://www.postgresql.org/docs/>`_ for installation instructions. Quick setup instructions, suitable for experimentation but not serious work, can be found in an `RDKit blog post <https://greglandrum.github.io/rdkit-blog/posts/2024-10-31-lwreg-and-the-cartridge.html>`_.



Configuration
-------------
If you plan to use the command line interface to lwreg, you need to set up a configuration file before you get started.
The configuration file is a JSON file that contains the following information:

    - **dbname**: The name of the database.
    - **dbtype**: The type of the database (sqlite or postgres).
    - **standardization**: The standardization options for the chemical structures.
    - **removeHs**: Whether to remove hydrogens from the chemical structures.
    - **useTautomerHashv2**: Whether to use the tautomer hash v2.
    - **registerConformers**: Whether to register conformers.
    - **numConformerDigits**: The number of digits to use when hashing conformer coordinates.
    - **lwregSchema**: The schema name for the lwreg tables (PostgreSQL only).

An example configuration file is shown below:

.. code-block:: json

    {
        "dbname": "test.db",
        "dbtype": "sqlite",
        "standardization": "fragment",
        "removeHs": 1,
        "useTautomerHashv2": 0,
        "registerConformers": 0
    }

If you are using lwreg through the Python API, you can pass the configuration as a dictionary to the lwreg functions.

Choosing the right standardization options for your project is crucial for the registration system to work properly.
There is a set of pre-defined standardization options including:

    - **none**: No standardization.
    - **sanitize**: Runs the standard RDKit sanitization on the molecule.
    - **fragment**: Generates the fragment parent of the molecule.
    - **charge**: Generates the charge parent of the molecule.
    - **tautomer**: Generates the tautomer parent of the molecule.
    - **super**: Generates the super parent of the molecule.
    - **canonicalize**: Canonicalizes the orientation of the molecule's 3D conformers (if present).

A user can also define their own standardization options. 

Besides the standardization options, there is also the possibility to define custom filers. For example, this filter rejects molecules with a net formal charge:

.. code-block:: python

    def reject_charged_molecules(mol):
        if Chem.GetFormalCharge(mol):
            return None
        return mol

Multiple standardization options and filters can be combined in a list in a user defined order.
The current collection of standardizers/filters is available in the module :code:`lwreg.standardization_lib`.
The chosen standardization pipeline is stored in the database itself. 

lwreg's Command Line interface
-------------------------------
Lwreg also provides a command line interface. ::

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


Running lwreg in Docker
-----------------------
lwreg can be run in a docker container. ::
    
    docker build -t lwreg .
    docker run -i -t -p 8888:8888 rdkit-lwreg /bin/bash -c "\
    apt update && apt install libtiff5 -y && \
    pip install notebook && \
    jupyter notebook \
    --notebook-dir=/lw-reg --ip='*' --port=8888 \
    --no-browser --allow-root"


Registering Conformers
----------------------
When the configuration option :code:`registerConformers` is set to True, lwreg expects that the compounds to be registered will have an associated conformer. 
The conformers are tracked in a different table than the molecule topologies and expectation is that every molecule registered will have a conformer (it's an error if they don't). 
It is possible to register multiple conformers for a single molecular structure (topology).
Note that once a database is created in :code:`registerConformers` mode, it probably should always be used in that mode. 
When in :code:`registerConformers` mode, the following behaviour in the API is changed:

- :code:`register()` and :code:`bulk_register()` require molecules to have associated conformers. Both return :code:`(molregno, conf_id)` tuples instead of just :code:`molregno` s.
- :code:`query()` can either be called with the :code:`ids` argument, which returns all of the conformers for the supplied molregnos as :code:`(molregno, conf_id)` tuples. If called with a molecule, the conformer of the molecule will be hashed and looked up in the conformers table, returning a list of :code:`(molregno,conf_id)` tuples.
- :code:`retrieve()` called with :code:`(molregno, conf_id)` tuples, it will return a dictionary of :code:`(molblock, 'mol')` tuples with :code:`(molregno, conf_id)` tuples as keys where the :code:`molblock`s contain the coordinates of the registered conformers.