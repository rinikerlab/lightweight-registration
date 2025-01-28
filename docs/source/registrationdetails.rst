Registration Details
=====================

The determination of uniqueness is the core of any chemical registration system, 
meaning the method to decide weather two molecular structures are the same or not.
Lwreg uses a hash-based approach to identify duplicates. 
As all comparisons are therefore reduced to string comparisons, lwreg is a fast approach in any modern database system.

Methodology of Hash Computation
-------------------------------
Lwreg makes uses of the RDKit's RegistrationHash functionality. 
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

The registration hash is then a SHA1 hash of all layers using Python's hashlib library.