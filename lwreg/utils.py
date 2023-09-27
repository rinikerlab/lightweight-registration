# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,

import rdkit
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import RegistrationHash
import csv
import json
import sqlite3
import enum
from . import standardization_lib

_violations = (sqlite3.IntegrityError, )
try:
    import psycopg2
    _violations = (sqlite3.IntegrityError, psycopg2.errors.UniqueViolation)
except ImportError:
    psycopg2 = None

from collections import namedtuple

standardizationOptions = {
    'none': lambda x: x,
    'sanitize': standardization_lib.RDKitSanitize(),
    'fragment': standardization_lib.FragmentParent(),
    'charge': standardization_lib.ChargeParent(),
    'tautomer': standardization_lib.TautomerParent(),
    'super': standardization_lib.SuperParent(),
}

_defaultConfig = json.loads('''{
    "dbname": "./testdb.sqlt",
    "dbtype": "sqlite3",
    "standardization": "fragment",
    "removeHs": 1,
    "use3DIfPresent": 1,
    "useTautomerHashv2": 0,
    "registerConformers": 0,
    "hashConformer": 0,
    "numConformerDigits": 3                            
}''')

from rdkit.Chem.RegistrationHash import HashLayer


class RegistrationFailureReasons(enum.Enum):
    DUPLICATE = enum.auto()
    PARSE_FAILURE = enum.auto()
    FILTERED = enum.auto()


def defaultConfig():
    return _defaultConfig.copy()


_config = {}


def _configure(filename='./config.json'):
    global _config
    if not _config:
        if filename:
            with open(filename, 'r') as inf:
                _config = json.load(inf)
        else:
            _config = defaultConfig()
    return _config


def _lookupWithDefault(config, key, defaults=_defaultConfig):
    return config.get(key, defaults[key])


_replace_placeholders_noop = lambda x: x
_replace_placeholders_pcts = lambda x: x.replace('?', '%s').replace('"', '')
_replace_placeholders = _replace_placeholders_noop


def _connect(config):
    global _replace_placeholders
    global _dbtype
    cn = config.get('connection', None)
    dbtype = _lookupWithDefault(config, 'dbtype').lower()
    if not cn:
        dbnm = config['dbname']
        if dbtype == 'sqlite3':
            uri = False
            if dbnm.startswith('file::'):
                uri = True
            cn = sqlite3.connect(dbnm, uri=uri)
        elif dbtype in ('postgres', 'postgresql'):
            dbtype = 'postgresql'
            if psycopg2 is None:
                raise ValueError("psycopg2 package not installed")
            cn = psycopg2.connect(database=dbnm,
                                  host=config.get('host', None),
                                  user=config.get('user', None),
                                  password=config.get('password', None))
    _dbtype = dbtype
    if dbtype == 'postgresql':
        _replace_placeholders = _replace_placeholders_pcts
    else:
        _replace_placeholders = _replace_placeholders_noop

    return cn


MolTuple = namedtuple('MolTuple', ('mol', 'datatype', 'rawdata'))


def _process_molblock(molblock, config):
    mol = Chem.MolFromMolBlock(molblock,
                               sanitize=False,
                               removeHs=_lookupWithDefault(config, 'removeHs'))
    if mol.GetConformer().Is3D() and _lookupWithDefault(
            config, 'use3DIfPresent'):
        Chem.AssignStereochemistryFrom3D(mol)
    return mol


def _parse_mol(mol=None, molfile=None, molblock=None, smiles=None, config={}):
    if mol is not None:
        datatype = 'pkl'
        raw = mol.ToBinary(propertyFlags=Chem.PropertyPickleOptions.AllProps)
    elif smiles is not None:
        spp = Chem.SmilesParserParams()
        spp.sanitize = False
        spp.removeHs = _lookupWithDefault(config, 'removeHs')
        mol = Chem.MolFromSmiles(smiles, spp)
        datatype = 'smiles'
        raw = smiles
    elif molblock is not None:
        mol = _process_molblock(molblock, config)
        datatype = 'mol'
        raw = molblock
    elif molfile is not None:
        with open(molfile, 'r') as inf:
            molblock = inf.read()
        mol = _process_molblock(molblock, config)
        datatype = 'mol'
        raw = molblock
    else:
        return MolTuple(None, None, None)

    return MolTuple(mol, datatype, raw)


def _get_standardization_list(config):
    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)
    sopts = _lookupWithDefault(config, 'standardization')

    res = []
    if type(sopts) not in (list, tuple):
        sopts = (sopts, )
    return tuple(sopts)


def _get_standardization_label(config):
    sopts = _get_standardization_list(config)
    res = []
    for sopt in sopts:
        if type(sopt) == str:
            nm = sopt
        elif hasattr(sopt, 'name'):
            nm = sopt.name
        else:
            nm = 'unknown'
        res.append(nm)
    return '|'.join(res)


def standardize_mol(mol, config=None):
    """ standardizes the input molecule using the 'standardization' 
    option in config and returns the result 

    The value of the 'standardization' option can be:
      - one of the strings from standardizationOptions
      - a callable object which should return either None (on failure) 
        or the standardized molecule (on success)
    A list or tuple of these is also ok, the entries in the list will 
    be applied in order.
    
    Keyword arguments:
    config -- configuration dict
    """
    sopts = _get_standardization_list(config)
    for sopt in sopts:
        if type(sopt) == str:
            sopt = standardizationOptions[sopt]
        mol = sopt(mol)
        if mol is None:
            return None
    return mol


def _get_conformer_hash(mol, numDigits):
    """ returns a hash for a molecule's conformer based on the atomic positions
    
    """
    if not mol.GetNumConformers():
        return ''
    ps = mol.GetConformer().GetPositions()
    hps = []
    for p in ps:
        coords = [str(round(x, numDigits)) for x in p]
        hps.append(','.join(coords))
    return ';'.join(sorted(hps))


def hash_mol(mol, escape=None, config=None):
    """ returns the hash layers and corresponding hash for the molecule 
    
    Keyword arguments:
    escape -- the escape layer
    config -- configuration dict
    """
    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)
    layers = RegistrationHash.GetMolLayers(
        mol,
        escape=escape,
        enable_tautomer_hash_v2=_lookupWithDefault(config,
                                                   'useTautomerHashv2'))

    if _lookupWithDefault(config, 'hashConformer'):
        escape = layers.get(RegistrationHash.HashLayer.ESCAPE, '')
        if escape:
            escape += '|\n'
        layers[
            RegistrationHash.HashLayer.ESCAPE] = escape + _get_conformer_hash(
                mol, _lookupWithDefault(config, "numConformerDigits"))

    mhash = RegistrationHash.GetMolHash(layers)

    return mhash, layers


def _register_mol(tpl,
                  escape,
                  cn,
                  curs,
                  config,
                  failOnDuplicate,
                  def_rdkit_version_label=None,
                  def_std_label=None):
    """ does the work of registering one molecule """
    registerConformers = _lookupWithDefault(config, "registerConformers")

    if registerConformers and tpl.mol is not None \
        and not tpl.mol.GetNumConformers():
        raise ValueError(
            "attempt to register a molecule without conformers when registerConformers is set"
        )

    standardization_label = _get_standardization_label(config)
    if def_std_label is None:
        curs.execute(
            "select value from registration_metadata where key='standardization'"
        )
        def_std_label = curs.fetchone()[0]
    if standardization_label == def_std_label:
        standardization_label = None

    if def_rdkit_version_label is None:
        curs.execute(
            "select value from registration_metadata where key='rdkitVersion'")
        def_rdkit_version_label = curs.fetchone()[0]
    rdkit_version_label = rdkit.__version__
    if rdkit_version_label == def_rdkit_version_label:
        rdkit_version_label = None
    mrn = None
    conf_id = None
    try:
        sMol = standardize_mol(tpl.mol, config=config)
        if sMol is None:
            return None, None
        molb = Chem.MolToV3KMolBlock(sMol)

        mhash, layers = hash_mol(sMol, escape=escape, config=config)

        regtuple = [mhash] + [
            layers[HashLayer.FORMULA], layers[HashLayer.CANONICAL_SMILES],
            layers[HashLayer.NO_STEREO_SMILES],
            layers[HashLayer.TAUTOMER_HASH],
            layers[HashLayer.NO_STEREO_TAUTOMER_HASH],
            layers[HashLayer.ESCAPE], layers[HashLayer.SGROUP_DATA]
        ]
        regtuple.append(rdkit_version_label)
        regtuple = tuple(regtuple)
        qs = ','.join('?' * len(regtuple))
        # will fail if the fullhash is already there
        if _dbtype != 'postgresql':
            curs.execute(
                _replace_placeholders(
                    f'insert into hashes values (NULL,{qs})'), regtuple)
            curs.execute('select last_insert_rowid()')
            mrn = curs.fetchone()[0]
        else:
            # Note that on postgresql the serial ids are increasing, but not sequential
            #  i.e. failed inserts will increment the counter
            curs.execute(
                _replace_placeholders(
                    f'insert into hashes values (default,{qs}) returning molregno'
                ), regtuple)
            mrn = curs.fetchone()[0]

        curs.execute(
            _replace_placeholders('insert into orig_data values (?, ?, ?)'),
            (mrn, tpl.rawdata, tpl.datatype))
        curs.execute(
            _replace_placeholders('insert into molblocks values (?, ?, ?)'),
            (mrn, molb, standardization_label))

        cn.commit()
    except _violations:
        cn.rollback()
        if failOnDuplicate and not (registerConformers
                                    and sMol.GetNumConformers()):
            raise
        else:
            curs.execute(
                _replace_placeholders(
                    'select molregno from hashes where fullhash=?'), (mhash, ))
            mrn = curs.fetchone()[0]
    except:
        cn.rollback()
        raise

    if registerConformers and sMol.GetNumConformers():
        try:
            chash = _get_conformer_hash(
                sMol, _lookupWithDefault(config, "numConformerDigits"))
            regtuple = (mrn, chash, molb)
            qs = '?,?,?'
            if _dbtype != 'postgresql':
                curs.execute(
                    _replace_placeholders(
                        f'insert into conformers values (NULL,{qs})'),
                    regtuple)
                curs.execute('select last_insert_rowid()')
                conf_id = curs.fetchone()[0]
            else:
                # Note that on postgresql the serial ids are increasing, but not sequential
                #  i.e. failed inserts will increment the counter
                curs.execute(
                    _replace_placeholders(
                        f'insert into conformers values (default,{qs}) returning conf_id'
                    ), regtuple)
                conf_id = curs.fetchone()[0]
            cn.commit()
        except _violations:
            cn.rollback()
            if failOnDuplicate:
                raise
            else:
                curs.execute(
                    _replace_placeholders(
                        'select conf_id from conformers where conformer_hash=?'
                    ), (chash, ))
                conf_id = curs.fetchone()[0]
        except:
            cn.rollback()
            raise

    return mrn, conf_id


def _get_delimiter(smilesfile):
    with open(smilesfile, 'r') as input_file:
        lines = input_file.readlines()
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(lines[-1])
        delimiter = dialect.delimiter
        expected_delimiters = ['\t', ',', ';', ' ', '  ', '   ', '    ', '|']
        if delimiter in expected_delimiters:
            return delimiter
        else:
            return


def _get_smiles_column(smilesfile, delimiter):
    with open(smilesfile, 'r', newline='') as input_file:
        rows = csv.reader(input_file, delimiter=delimiter)
        rows = list(rows)
        for row in rows[::-1]:
            with rdBase.BlockLogs() as bl:
                for i, smi in enumerate(row):
                    if Chem.MolFromSmiles(smi, sanitize=False) is not None:
                        return i
        raise ValueError('No valid SMILES found in file.')


def _has_header(smilesfile, delimiter):
    with open(smilesfile, 'r') as input_file:
        header_line = input_file.readline()[:-1]
        elements = header_line.split(delimiter)
        with rdBase.BlockLogs() as bl:
            for element in elements:
                if Chem.MolFromSmiles(element, sanitize=False):
                    return False
        return True


def _get_mols_from_smilesfile(smilesfile):
    delimiter = _get_delimiter(smilesfile)
    if not delimiter:
        col = 0
        delimiter = ' '
    else:
        col = _get_smiles_column(smilesfile, delimiter)
    has_header = _has_header(smilesfile=smilesfile, delimiter=delimiter)
    if not has_header:
        with open(smilesfile, 'r') as input_file:
            content = 'DUMMY_HEADER\n' + input_file.read()
        mols = Chem.SmilesMolSupplierFromText(content,
                                              delimiter=delimiter,
                                              smilesColumn=col)
    else:
        mols = Chem.SmilesMolSupplier(smilesfile,
                                      delimiter=delimiter,
                                      smilesColumn=col)
    return mols


def register(config=None,
             mol=None,
             molfile=None,
             molblock=None,
             smiles=None,
             escape=None,
             fail_on_duplicate=True,
             no_verbose=True):
    """ registers a new molecule, assuming it doesn't already exist,
    and returns the new registry number (molregno)
    
    only one of the molecule format objects should be provided

    Keyword arguments:
    config     -- configuration dict
    mol        -- RDKit molecule object
    molfile    -- MOL or SDF filename
    molblock   -- MOL or SDF block
    smiles     -- smiles
    escape     -- the escape layer
    failOnDuplicate -- if true then an exception is raised when trying to register a duplicate
    no_verbose -- if this is False then the registry number will be printed
    """
    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)
    tpl = _parse_mol(mol=mol,
                     molfile=molfile,
                     molblock=molblock,
                     smiles=smiles,
                     config=config)

    cn = _connect(config)
    curs = cn.cursor()
    if tpl.mol is None:
        return RegistrationFailureReasons.PARSE_FAILURE

    mrn, conf_id = _register_mol(tpl, escape, cn, curs, config,
                                 fail_on_duplicate)
    if mrn is None:
        return RegistrationFailureReasons.FILTERED
    if not _lookupWithDefault(config, "registerConformers"):
        res = mrn
    else:
        res = mrn, conf_id
    if not no_verbose:
        print(res)
    return res


def bulk_register(config=None,
                  mols=None,
                  sdfile=None,
                  smilesfile=None,
                  escapeProperty=None,
                  failOnDuplicate=True,
                  no_verbose=True):
    """ registers multiple new molecules, assuming they don't already exist,
    and returns the new registry numbers (molregno)
    
    The result tuple includes a single entry for each molecule in the input.
    That entry can be one of the following:
      - the registry number (molregno) of the registered molecule
      - RegistrationFailureReasons.DUPLICATE if failOnDuplicate is True and a
        molecule is a duplicate
      - RegistrationFalureReasons.PARSE_FAILURE if there was a problem processing
        the molecule 
    
    only one of the molecule format objects should be provided

    Keyword arguments:
    config         -- configuration dict
    mols           -- an iterable of RDKit molecule objects
    sdfile         -- SDF filename
    escapeProperty -- the molecule property to use as the escape layer
    failOnDuplicate -- if true then RegistraionFailureReasons.DUPLICATE will be returned 
                       for each already-registered molecule, otherwise the already existing
                       structure ID will be returned
    no_verbose     -- if this is False then the registry numbers will be printed
    """

    if mols:
        pass
    elif sdfile:
        mols = Chem.ForwardSDMolSupplier(sdfile)
    elif smilesfile:
        mols = _get_mols_from_smilesfile(smilesfile)
    else:
        raise ValueError('No input molecules provided!')

    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)
    res = []
    cn = _connect(config)
    curs = cn.cursor()

    curs.execute(
        "select value from registration_metadata where key='standardization'")
    def_std_label = curs.fetchone()[0]
    curs.execute(
        "select value from registration_metadata where key='rdkitVersion'")
    def_rdkit_version_label = curs.fetchone()[0]

    for mol in mols:
        if mol is None:
            res.append(RegistrationFailureReasons.PARSE_FAILURE)
            continue
        tpl = _parse_mol(mol=mol, config=config)
        try:
            if escapeProperty is not None and mol.HasProp(escapeProperty):
                escape = mol.GetProp(escapeProperty)
            else:
                escape = None
            mrn, conf_id = _register_mol(
                tpl,
                escape,
                cn,
                curs,
                config,
                failOnDuplicate,
                def_rdkit_version_label=def_rdkit_version_label,
                def_std_label=def_std_label)

            if mrn is None:
                mrn = RegistrationFailureReasons.FILTERED

            if not _lookupWithDefault(config, "registerConformers"):
                res.append(mrn)
            else:
                res.append((mrn, conf_id))
        except _violations:
            res.append(RegistrationFailureReasons.DUPLICATE)
    if not no_verbose:
        print(res)
    return tuple(res)


def _getConfIdsForMolregnos(ids, config):
    assert _lookupWithDefault(config, "registerConformers")

    cn = _connect(config)
    curs = cn.cursor()

    qs = _replace_placeholders(','.join('?' * len(ids)))
    curs.execute(
        f'select molregno,conf_id from conformers where molregno in ({qs})',
        ids)
    res = curs.fetchall()
    return res


def query(config=None,
          layers='ALL',
          ids=None,
          mol=None,
          molfile=None,
          molblock=None,
          smiles=None,
          escape=None,
          no_verbose=True):
    """ queries to see if a molecule has already been registered,
    and returns the corresponding registry numbers (molregnos)
    
    only one of the molecule format objects should be provided

    Keyword arguments:
    config     -- configuration dict
    layers     -- hash layers to be used to determine identity
    ids        -- list or tuple of molregnos. 
                  Only makes sense if registerConformers is set, 
                  in which case this will return all conf_ids for 
                  the molregnos in the ids list as a list of 
                  (molregno, conf_id) tuples
    mol        -- RDKit molecule object
    molfile    -- MOL or SDF filename
    molblock   -- MOL or SDF block
    smiles     -- smiles
    escape     -- the escape layer
    no_verbose -- if this is False then the registry numbers will be printed
    """
    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)

    if ids is not None:
        if not _lookupWithDefault(config, "registerConformers"):
            raise ValueError(
                'passing ids to query() only makes sense when registerConformers is enabled'
            )
        res = _getConfIdsForMolregnos(ids, config=config)
    else:
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
            curs.execute(
                _replace_placeholders(
                    'select molregno from hashes where fullhash=?'), (mhash, ))
        else:
            vals = []
            query = []
            if type(layers) == str:
                layers = layers.upper().split(',')
            if not hasattr(layers, '__len__'):
                layers = [layers]
            for lyr in layers:
                if type(lyr) == str:
                    k = getattr(HashLayer, lyr)
                else:
                    k = lyr
                    lyr = str(lyr).split('.')[-1]
                vals.append(hlayers[k])
                query.append(f'"{lyr}"=?')

            query = _replace_placeholders(' and '.join(query))
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
    """ returns the molecule data for one or more registry ids (molregnos)
    The return value is a tuple of (molregno, data, format) 3-tuples    


    only one of id or ids should be provided

    Keyword arguments:
    config       -- configuration dict
    ids          -- an iterable of registry ids (molregnos)
    id           -- a registry id (molregno)
    as_submitted -- if True, then the structure will be returned as registered
    no_verbose   -- if this is False then the registry number will be printed
    """
    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)
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
    qs = _replace_placeholders(','.join('?' * len(ids)))
    curs.execute(f'select {qry} where molregno in ({qs})', ids)

    res = curs.fetchall()
    if not no_verbose:
        if res:
            for entry in res:
                print(entry)
        else:
            print('not found')

    return tuple(res)


def _registerMetadata(curs, config):
    dc = defaultConfig()
    dc.update(config)
    for k, v in dc.items():
        curs.execute(
            _replace_placeholders(
                'insert into registration_metadata values (?,?)'),
            (str(k), str(v)))

    curs.execute(
        _replace_placeholders(
            'insert into registration_metadata values (?,?)'),
        ('rdkitVersion', rdkit.__version__))


def initdb(config=None, confirm=False):
    """ initializes the registration database    

    NOTE that this call destroys any existing information in the registration database

    Keyword arguments:
    config  -- configuration dict
    confirm -- if set to False we immediately return
    """
    if not confirm:
        return
    if not config:
        config = _configure()
    elif type(config) == str:
        config = _configure(filename=config)
    cn = _connect(config)
    curs = cn.cursor()

    curs.execute('drop table if exists registration_metadata')
    curs.execute('create table registration_metadata (key text, value text)')
    _registerMetadata(curs, config)

    curs.execute('drop table if exists hashes')
    if _dbtype != 'postgresql':
        curs.execute(
            '''create table hashes (molregno integer primary key, fullhash text unique, 
            formula text, canonical_smiles text, no_stereo_smiles text, 
            tautomer_hash text, no_stereo_tautomer_hash text, "escape" text, sgroup_data text, rdkitVersion text)'''
        )
    else:
        curs.execute(
            '''create table hashes (molregno serial primary key, fullhash text unique, 
            formula text, canonical_smiles text, no_stereo_smiles text, 
            tautomer_hash text, no_stereo_tautomer_hash text, "escape" text, sgroup_data text, rdkitVersion text)'''
        )
    curs.execute('drop table if exists orig_data')
    curs.execute(
        'create table orig_data (molregno integer primary key, data text, datatype text)'
    )
    curs.execute('drop table if exists molblocks')
    curs.execute(
        'create table molblocks (molregno integer primary key, molblock text, standardization text)'
    )

    curs.execute('drop table if exists conformers')
    if _lookupWithDefault(config, "registerConformers"):
        if _dbtype != 'postgresql':
            curs.execute(
                '''create table conformers (conf_id integer primary key, molregno integer not null, 
                   conformer_hash text not null unique, molblock text)''')
        else:
            curs.execute(
                '''create table conformers (conf_id serial primary key, molregno integer not null, 
                   conformer_hash text not null unique, molblock text)''')

    cn.commit()
    return True
