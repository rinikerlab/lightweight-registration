#! /bin/env python
import click
import sqlite3
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import RegistrationHash

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
def initdb():
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
        'create table hashes (molregno int primary key, fullhash text unique, formula text, smiles text, nochi_smiles text, tautomer text, nochi_tautomer text, escape_layer text, sgroup_data text)'
    )


@cli.command()
@click.option(
    "--smiles",
    default=None,
)
@click.option(
    "--escape",
    default=None,
)
def register(molfile=None, molblock=None, smiles=None, escape=None):
    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        datatype = 'smiles'
        ind = smiles
        molb = Chem.MolToV3KMolBlock(mol)
    cn = _connect()
    mrn = _getNextRegno(cn)
    curs = cn.cursor()
    curs.execute('insert into orig_data values (?, ?, ?)',
                 (mrn, ind, datatype))
    curs.execute('insert into molblocks values (?, ?)', (mrn, molb))
    sMol = rdMolStandardize.FragmentParent(mol)
    layers = RegistrationHash.GetMolLayers(sMol,escape=escape)
    mhash = RegistrationHash.GetMolHash(layers)

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


@cli.command()
@click.option('--who', default='world')
def greet(who):
    print(f'hello {who}')


if __name__ == '__main__':
    cli()