Quick Start
===========

.. _GetStarted:


What is a registration system?
------------------------------
A chemical registration system is a software system (e.g. database) which holds chemical structures and assigns them identifiers.
In order to do so, a chemical registration system needs to be able to check a presented compound against already registered ones and decide upon its uniqueness. 
The identifiers can be used as keys when storing additional data about the chemical structures in additional tables.

Lwreg's Basic Functionality
---------------------------
There are four basic operations provides by the lwreg package that can either by called via the Python API or the command line interface:

.. list-table:: Basic operations
   :widths: 10 30
   :header-rows: 1

   * - Operation
     - Description
   * - register
     - Register a new compound, returns a unique identifier (molregno) upon success.
   * - query
     - Query the database if a compound is already registered.
       If so, returns the respective molregno(s) matching the search.
   * - retrieve
     - Get the resgisted structure(s) for one or multiple given identifier(s).

See the :ref:`Basic Functionality Reference` for more details.

Installation Guide
------------------
The following outline assumes that the user has conda (or an equivalent) installed. 

.. code-block:: bash

    conda env create --name py311_lwreg --file=https://raw.githubusercontent.com/rinikerlab/lightweight-registration/main/environment.yml
    conda activate py311_lwreg
    python -m pip install git+https://github.com/rinikerlab/lightweight-registration

To verify the installation by running

.. code-block:: bash

    lwreg --help

When starting to use lwreg, you also must choose your prefered database management system. 
The two systems compatible with lwreg are SQLite and PostgreSQL.

To check if SQLite is available in your Python environment, run

.. code-block:: bash

    python -c "import sqlite3; print(sqlite3.sqlite_version)"

For setting up a PostgreSQL database, you need to have the PostgreSQL server installed and running.
Please refer to the `PostgreSQL documentation <https://www.postgresql.org/docs/>`_ for installation instructions.
If you choose to use PostgreSQL, you also need to install the psycopg2 package:

.. code-block:: bash

    conda install -c conda-forge psycopg2

Configuration
-------------
Before the lwreg initialization, you need to set up a configuration file.
The configuration file is a YAML file that containing the following information:
The configuration file is a YAML file that contains the following information:

    - **dbname**: The name of the database.
    - **dbtype**: The type of the database (sqlite or postgres).
    - **standardization**: The standardization options for the chemical structures.
    - **removeHs**: Whether to remove hydrogens from the chemical structures.
    - **useTautomerHashv2**: Whether to use the tautomer hash v2.
    - **registerConformers**: Whether to register conformers.
    - **numConformerDigits**: The number of digits to use when hashing conformer coordinates.
    - **lwregSchema**: The schema name for the lwreg tables (PostgreSQL only).

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

Besides the standardization options, there is also the possibility to define custom filers:

.. code-block:: python

    def reject_charged_molecules(mol):
        if Chem.GetFormalCharge(mol):
            return None
        return mol

Multiple standardization options and filters can be combined in a list in a user defined order.
The chosen standardization pipeline is stored in the database itself. 

Registering Conformers
----------------------
