#! /bin/env python

# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE, 
import click
import sqlite3
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import RegistrationHash
import utils
import json
_config = {}


def _connect():
    cn = sqlite3.connect(_config['dbfile'])
    return cn


def _configure(filename='./config.json'):
    global _config
    with open(filename, 'r') as inf:
        _config = json.load(inf)


def _getNextRegno(cn):
    curs = cn.cursor()
    curs.execute('select max(molregno) from orig_data')
    row = curs.fetchone()
    if row[0] is None:
        res = 1
    else:
        res = row[0] + 1
    return res


@click.group()
def cli():
    _configure()


@cli.command()
@click.option(
    "--confirm",
    default='no',
)
def initdb(confirm='no'):
    if confirm != 'yes':
        click.echo("inidb not confirmed, aborting")
        return

    print(f'dbname: {_config["dbfile"]}')
    cn = _connect()
    curs = cn.cursor()
    curs.execute('drop table if exists orig_data')
    curs.execute(
        'create table orig_data (molregno int primary key, data text, datatype text)'
    )
    curs.execute('drop table if exists molblocks')
    curs.execute(
        'create table molblocks (molregno int primary key, molblock text)')
    curs.execute('drop table if exists hashes')
    curs.execute(
        '''create table hashes (molregno int primary key, fullhash text unique, 
          formula text, canonical_smiles text, no_stereo_smiles text, 
          tautomer_hash text, no_stereo_tautomer_hash text, "escape" text, sgroup_data text)'''
    )



@cli.command()
@click.option(
    "--ids",
    default=None,
)
@click.option(
    "--id",
    default=None,
)
@click.option(
    "--as_submitted",
    default=False,
    is_flag=True
)
@click.option(
    "--no-verbose",
    default=False,
    is_flag=True
)
def retrieve(ids=None, id=None, as_submitted=False, no_verbose=True):
    if id is not None:
        ids = [int(id)]
    if ids is not None:
        if type(ids)==str:
            ids = [int(x) for x in ids.split(',')]
    cn = _connect()
    curs = cn.cursor()
    if as_submitted:
        qry = 'molregno,data,datatype from orig_data'
    else:
        qry = "molregno,molblock,'mol' from molblocks"
    qs = ','.join('?'*len(ids))
    curs.execute(f'select {qry} where molregno in ({qs})',ids)
    
    res = curs.fetchall()
    if not no_verbose:
        if res:
            for entry in res:
                print(entry)
        else:
            print('not found')

    return res


@cli.command()
@click.option(
    "--smiles",
    default=None,
)
@click.option(
    "--escape",
    default=None,
)
@click.option(
    "--layers",
    default='ALL',
)
@click.option(
    "--no-verbose",
    default=False,
    is_flag=True
)
def query(layers='ALL',molfile=None, molblock=None, smiles=None, escape=None, no_verbose=True):
    tpl = utils.parse_mol(molfile=molfile,molblock=molblock,smiles=smiles,config=_config)
    sMol = utils.standardize_mol(tpl.mol,config=_config)
    mhash,hlayers = utils.hash_mol(sMol,escape=escape,config=_config)

    cn = _connect()
    curs = cn.cursor()
    layers = layers.upper()
    if layers=='ALL':
        curs.execute('select molregno from hashes where fullhash=?',(mhash,))
    else:
        vals = []
        query = []
        if type(layers)==str:
            layers = layers.split(',')
        for lyr in layers:
            k = getattr(RegistrationHash.HashLayer,lyr)
            vals.append(hlayers[k])
            query.append(f'"{lyr}"=?')

        query = ' and '.join(query)
        curs.execute(f'select molregno from hashes where {query}',vals)
    
    res = [x[0] for x in curs.fetchall()]
    if not no_verbose:
        if res:
            print(' '.join(str(x) for x in res))
        else:
            print('not found')

    return res

@cli.command()
@click.option(
    "--smiles",
    default=None,
)
@click.option(
    "--escape",
    default=None,
)
@click.option(
    "--no-verbose",
    default=False,
    is_flag=True
)
def register(molfile=None, molblock=None, smiles=None, escape=None, no_verbose=True):
    tpl = utils.parse_mol(molfile=molfile,molblock=molblock,smiles=smiles,config=_config)
    
    molb = Chem.MolToV3KMolBlock(tpl.mol)
    cn = _connect()
    mrn = _getNextRegno(cn)
    curs = cn.cursor()
    curs.execute('insert into orig_data values (?, ?, ?)',
                 (mrn, tpl.rawdata, tpl.datatype))
    curs.execute('insert into molblocks values (?, ?)', (mrn, molb))

    sMol = utils.standardize_mol(tpl.mol,config=_config)

    mhash,layers = utils.hash_mol(sMol,escape=escape,config=_config)

    # will fail if the fullhash is already there
    curs.execute('insert into hashes values (?,?,?,?,?,?,?,?,?)', (
        mrn,
        mhash,
        layers[RegistrationHash.HashLayer.FORMULA],
        layers[RegistrationHash.HashLayer.CANONICAL_SMILES],
        layers[RegistrationHash.HashLayer.NO_STEREO_SMILES],
        layers[RegistrationHash.HashLayer.TAUTOMER_HASH],
        layers[RegistrationHash.HashLayer.NO_STEREO_TAUTOMER_HASH],
        layers[RegistrationHash.HashLayer.ESCAPE],
        layers[RegistrationHash.HashLayer.SGROUP_DATA],
    ))

    cn.commit()
    if not no_verbose:
        print(mrn)
    return mrn

@cli.command()
@click.option('--who', default='world')
def greet(who):
    print(f'hello {who}')


if __name__ == '__main__':
    cli()