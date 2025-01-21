# Copyright (C) 2022-2025 Greg Landrum and other lwreg contributors
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,

from . import utils

import logging
try:
    import psycopg2
except ImportError:
    psycopg2 = None
    logging.INFO("psycopg2 not available, RDKit cartridge support disabled")


# decorator to disable functions if psycopg2 is not available
def psycopg2_available(func):

    def pass_through(*args, **kwargs):
        if psycopg2 is None:
            raise NotImplementedError("psycopg2 is not available")
        return func(*args, **kwargs)

    return pass_through


@psycopg2_available
def enable_cartridge(config):
    conn = utils.connect(config)
    curs = conn.cursor()
    curs.execute("create extension if not exists rdkit")


@psycopg2_available
def disable_cartridge(config):
    conn = utils.connect(config)
    curs = conn.cursor()
    curs.execute("drop extension if exists rdkit")


@psycopg2_available
def populate_rdkit_schema(config, force=False):
    if not force:
        print("This will destroy any existing information in the rdkit schema")
        response = input("  are you sure? [yes/no]: ")
        if response != 'yes':
            print("cancelled")
            return False

    enable_cartridge(config)
    conn = utils.connect(config)
    curs = conn.cursor()
    rdkit_schema_name = config.get('rdkit_schema', 'rdk')
    lwreg_schema_name = config.get('lwregSchema', 'public')
    if not lwreg_schema_name:
        lwreg_schema_name = 'public'
    curs.execute(f"create schema if not exists {rdkit_schema_name}")
    curs.execute(f"drop table if exists {rdkit_schema_name}.mols cascade")
    curs.execute(
        f"create table {rdkit_schema_name}.mols (molregno integer not null unique references {lwreg_schema_name}.hashes (molregno), m mol)"
    )
    curs.execute(
        f"insert into {rdkit_schema_name}.mols select molregno, mol_from_ctab(molblock::cstring,false) m from {lwreg_schema_name}.molblocks"
    )
    curs.execute(
        f'create index {rdkit_schema_name}_molidx on {rdkit_schema_name}.mols using gist(m)'
    )
    curs.execute(
        f'''CREATE OR REPLACE FUNCTION {rdkit_schema_name}_copy_new_mol() RETURNS TRIGGER AS
$BODY$
BEGIN
    INSERT INTO
        {rdkit_schema_name}.mols(molregno,m)
        VALUES(new.molregno,mol_from_ctab(new.molblock::cstring,false));
           RETURN new;
END;
$BODY$
language plpgsql;''')
    curs.execute(
        f'''CREATE OR REPLACE TRIGGER {lwreg_schema_name}_mol_insert_trigger
     AFTER INSERT ON {lwreg_schema_name}.molblocks
     FOR EACH ROW
     EXECUTE FUNCTION {rdkit_schema_name}_copy_new_mol();''')
    conn.commit()
    return True
