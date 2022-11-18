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