# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import RegistrationHash
import json
import sqlite3

from collections import namedtuple

_config = {}


def _configure(filename='./config.json'):
    global _config
    if not _config:
        with open(filename, 'r') as inf:
            _config = json.load(inf)
    return _config


def _connect(config):
    cn = config.get('connection', None)
    if not cn:
        uri = False
        dbnm = config['dbfile']
        if dbnm.startswith('file::'):
            uri = True
        cn = sqlite3.connect(dbnm, uri=uri)
    return cn


def _getNextRegno(cn):
    curs = cn.cursor()
    curs.execute('select max(molregno) from orig_data')
    row = curs.fetchone()
    if row[0] is None:
        res = 1
    else:
        res = row[0] + 1
    return res


MolTuple = namedtuple('MolTuple', ('mol', 'datatype', 'rawdata'))


def _parse_mol(mol=None, molfile=None, molblock=None, smiles=None, config={}):
    if mol is not None:
        datatype = 'pkl'
        raw = mol.ToBinary()
    elif smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        datatype = 'smiles'
        raw = smiles
    elif molblock is not None:
        mol = Chem.MolFromMolBlock(molblock)
        datatype = 'mol'
        raw = molblock
    elif molfile is not None:
        with open(molfile, 'r') as inf:
            molblock = inf.read()
        mol = Chem.MolFromMolBlock(molblock)
        datatype = 'mol'
        raw = molblock

    return MolTuple(mol, datatype, raw)


def standardize_mol(mol, config=None):
    if config is None:
        config = _configure()
    sMol = rdMolStandardize.FragmentParent(mol)
    return sMol


def hash_mol(mol, escape=None, config=None):
    if config is None:
        config = _configure()
    layers = RegistrationHash.GetMolLayers(mol, escape=escape)
    mhash = RegistrationHash.GetMolHash(layers)
    return mhash, layers

def _register_mol(tpl,escape,cn,curs,config):
    molb = Chem.MolToV3KMolBlock(tpl.mol)
    mrn = _getNextRegno(cn)
    try:
        curs.execute('insert into orig_data values (?, ?, ?)',
                     (mrn, tpl.rawdata, tpl.datatype))
        curs.execute('insert into molblocks values (?, ?)', (mrn, molb))

        sMol = standardize_mol(tpl.mol, config=config)

        mhash, layers = hash_mol(sMol, escape=escape, config=config)

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
    except:
        cn.rollback()
        raise
    return mrn

def register(config=None,
             mol=None,
             molfile=None,
             molblock=None,
             smiles=None,
             escape=None,
             no_verbose=True):
    if config is None:
        config = _configure()
    tpl = _parse_mol(mol=mol,
                     molfile=molfile,
                     molblock=molblock,
                     smiles=smiles,
                     config=config)

    cn = _connect(config)
    curs = cn.cursor()
    mrn = _register_mol(tpl,escape,cn,curs,config)
    if not no_verbose:
        print(mrn)
    return mrn


def query(config=None,
          layers='ALL',
          mol=None,
          molfile=None,
          molblock=None,
          smiles=None,
          escape=None,
          no_verbose=True):
    if config is None:
        config = _configure()

    tpl = _parse_mol(mol=mol,
                     molfile=molfile,
                     molblock=molblock,
                     smiles=smiles,
                     config=config)
    sMol = standardize_mol(tpl.mol, config=config)
    mhash, hlayers = hash_mol(sMol, escape=escape, config=config)

    cn = _connect(config)
    curs = cn.cursor()
    if layers == 'ALL':
        curs.execute('select molregno from hashes where fullhash=?', (mhash, ))
    else:
        vals = []
        query = []
        if type(layers) == str:
            layers = layers.upper().split(',')
        for lyr in layers:
            if type(lyr) == str:
                k = getattr(RegistrationHash.HashLayer, lyr)
            else:
                k = lyr
                lyr = str(lyr).split('.')[-1]
            vals.append(hlayers[k])
            query.append(f'"{lyr}"=?')

        query = ' and '.join(query)
        curs.execute(f'select molregno from hashes where {query}', vals)

    res = [x[0] for x in curs.fetchall()]
    if not no_verbose:
        if res:
            print(' '.join(str(x) for x in res))
        else:
            print('not found')

    return res


def retrieve(config=None,
             ids=None,
             id=None,
             as_submitted=False,
             no_verbose=True):
    if config is None:
        config = _configure()
    if id is not None:
        ids = [int(id)]
    if ids is not None:
        if type(ids) == str:
            ids = [int(x) for x in ids.split(',')]
    cn = _connect(config)
    curs = cn.cursor()
    if as_submitted:
        qry = 'molregno,data,datatype from orig_data'
    else:
        qry = "molregno,molblock,'mol' from molblocks"
    qs = ','.join('?' * len(ids))
    curs.execute(f'select {qry} where molregno in ({qs})', ids)

    res = curs.fetchall()
    if not no_verbose:
        if res:
            for entry in res:
                print(entry)
        else:
            print('not found')

    return res


def initdb(config=None):
    if config is None:
        config = _configure()
    cn = _connect(config)
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
    return True