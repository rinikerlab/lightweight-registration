Registration Details
=====================

The method for the determination of uniqueness, deciding whether or not two structures are the same, is the core of any chemical registration system. lwreg uses a hash-based approach to identify duplicates. 
As all comparisons are simply string comparisons, lwreg is a fast approach in any modern database system and does not require any database extensions for handling chemical structures.

Molecule Hash Computation
--------------------------
lwreg makes uses of the `RDKit's RegistrationHash <https://rdkit.org/docs/source/rdkit.Chem.RegistrationHash.html>`_ functionality.
It consists of seven hash layers:


.. list-table:: Hash layers
   :widths: 3 10 30
   :header-rows: 1

   * - Hash layer number
     - Layer name
     - Description
   * - 1
     - FORMULA 
     - molecular formula
   * - 2
     - CANONICAL_SMILES
     - canonical SMILES representation
   * - 3
     - TAUTOMER_HASH
     - HetAtomTautomer hash (either v1 or v2) from MolHash. This includes information about stereochemistry (including enhanced stereochemistry)
   * - 4
     - NO_STEREO_SMILES
     - canonical CXSMILES for the molecule without any information about stereochemistry
   * - 5
     - NO_STEREO_TAUTOMER_HASH
     - HetAtomTautomer hash (either v1 or v2) for the molecule without any information about stereochemistry
   * - 6
     - SGROUP_DATA
     - canonicalized form of some of the molecule's SGroup data (if any). Which SGroup data fields are used can be configured by the user
   * - 7
     - ESCAPE
     - free-text field allowing arbitrary information to be used as part of the molecule hash

The registration hash is then a SHA1 hash of all layers computed using Python's hashlib library.

Conformer Hash Computation
---------------------------

We use a simple hashing scheme to quickly recognize whether or not a particular conformer has already been seen in the database. The algorithm used for hashing a conformer is:

* Convert the (x,y,z) coordinates of each atom into strings by rounding each coordinate to a particular number of digits after the decimal (default: 3).
* Describe the position of each atom as a string by concatenating the strings of the individual coordinates together with commas.
* Sort the string representations of the atom positions lexicographically.
* Construct the final hash by concatenating the sorted string representations together with semicolons and generating the SHA256 hash of the result.

Note that this simple scheme is independent of the atom ordering but it is neither translationally nor rotationally invariant. This is essential to allow the system to be used for storing pre-aligned conformers (e.g., for storing docking poses). If the user desires translational and/or rotational invariance, they should either standardize the orientation of conformers themselves before registering them or use the ``CanonicalizeOrientation`` step in the standardization pipeline to automate that process.

If registering a multi-conformer molecule, it is most efficient to call :code:`register_multiple_conformers()`. 
That only does the work of standardizing the molecule and calculating the molecule hash once.

The Base Tables
----------------

Lwreg will create the following tables in the database:

    - **registration_metadata**: Stores metadata about the registration system.
    - **hashes**: Stores the hashes for each molecule. 
        - Fields: molregno, fullhash, formula, canonical_smiles, no_stereo_smiles, tautomer_hash, no_stereo_tautomer_hash, escape, sgroup_data, rdkitVersion
    - **orig_data**: Stores the original data for each molecule.
        - Fields: molregno, data, datatype, timestamp
    - **molblocks**: Stores the molblocks for each molecule.
        - Fields: molregno, molblock, standardization
    - **conformers (only created in conformer mode)**: Stores the conformers for each molecule. 
        - Fields: conf_id, molregno, conformer_hash, molblock