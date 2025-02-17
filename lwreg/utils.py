# Copyright (C) 2022-2024 Greg Landrum and other lwreg contributors
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,

import rdkit
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import RegistrationHash
from hashlib import sha256
import csv
import json
import sqlite3
import enum
import os.path
from . import standardization_lib
import logging
from tqdm import tqdm
import base64

_violations = (sqlite3.IntegrityError, )
try:
    import psycopg2
    _violations = (sqlite3.IntegrityError, psycopg2.errors.UniqueViolation)
except ImportError:
    psycopg2 = None

from collections import namedtuple

standardizationOptions = {
    'none': standardization_lib.NoStandardization(),
    'sanitize': standardization_lib.RDKitSanitize(),
    'fragment': standardization_lib.FragmentParent(),
    'charge': standardization_lib.ChargeParent(),
    'tautomer': standardization_lib.TautomerParent(),
    'super': standardization_lib.SuperParent(),
    'canonicalize': standardization_lib.CanonicalizeOrientation(),
}

_defaultConfig = {
    "dbname": "./testdb.sqlt",
    "dbtype": "sqlite3",
    "standardization": "fragment",
    "removeHs": 1,
    "useTautomerHashv2": 0,
    "registerConformers":
    0,  # toggle registering conformers as well as compound structures
    "numConformerDigits":
    3,  # number of digits to use when hashing conformer coordinates
    "lwregSchema":
    "",  # the schema name to use for the lwreg tables (no effect with sqlite3)
}

from rdkit.Chem.RegistrationHash import HashLayer
from rdkit.Chem import KekulizeException


class RegistrationFailureReasons(enum.Enum):
    DUPLICATE = enum.auto()
    PARSE_FAILURE = enum.auto()
    FILTERED = enum.auto()


def defaultConfig():
    return _defaultConfig.copy()


_config = {}


def configure_from_database(dbname=None,
                            connection=None,
                            dbtype=None,
                            host=None,
                            user=None,
                            password=None,
                            lwregSchema=None):
    ''' returns a config dict with values from the registration 
    metadata table in the database.

    This will overwrite some of the values already in the config dict.

    Note that in order for this to work the arguments must provide whatever
    information is needed to connect to the database. This can be 'connection'
    with a direct connection object or 'dbname' and 'dbtype' (potentially with
    'host', 'user', and 'password' if those are required). If you used a
    nondefault schema when initializing the database, you'll also need to
    provide 'lwregSchema' here.

    If 'dbtype' is not provided, the following heuristics are used:
      - if 'dbname' corresponds to an existing file, then sqlite3 is used
      - if 'host' is provided, then postgresql is used
      - otherwise the default dbtype, currently sqlite3, is used

    Keyword arguments: 
    dbname -- the name of the database (one of dbname or connection must be provided)
    connection -- a connection object (one of dbname or connection must be provided)
    dbtype -- the type of database (sqlite3 or postgresql)
    host -- the host to connect to (for postgresql)
    user -- the user to connect as (for postgresql)
    password -- the password to use (for postgresql)
    lwregSchema -- the schema name to use for the lwreg tables (for postgresql)
        
    '''
    global _config
    config = {}
    if connection is not None:
        config['connection'] = connection
    elif dbname is not None:
        config['dbname'] = dbname
    else:
        raise ValueError(
            'at least one of dbname or connection must be provided')

    if dbtype is not None:
        config['dbtype'] = dbtype
    else:
        if os.path.exists(dbname):
            # if the db is a file, then we'll assume sqlite
            config['dbtype'] = 'sqlite3'
        elif host is not None:
            # if they provided a host, it's probably postgresql
            config['dbtype'] = 'postgresql'

    if host is not None:
        config['host'] = host
    if user is not None:
        config['user'] = user
    if password is not None:
        config['password'] = password
    if lwregSchema is not None:
        config['lwregSchema'] = lwregSchema

    if 'connection' in config:
        cn = config['connection']
    else:
        cn = connect(config)
    curs = cn.cursor()
    curs.execute(f'select * from {registrationMetadataTableName}')
    rows = curs.fetchall()
    for k, v in rows:
        if k in ('rdkitVersion', 'user', 'password'):
            continue
        try:
            v = int(v)
        except ValueError:
            try:
                v = float(v)
            except ValueError:
                pass
        config[k] = v
    _config = config
    return config


def set_default_config(config):
    ''' sets the default configuration to be used by the other functions in 
    this module to the configuration object which is passed in
    
    Arguments:
    config -- configuration dict

    '''
    global _config
    assert isinstance(config, dict)
    _config = config


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

_baseregistrationMetadataTableName = 'registration_metadata'
_baseorigDataTableName = 'orig_data'
_basehashTableName = 'hashes'
_basemolblocksTableName = 'molblocks'
_baseconformersTableName = 'conformers'

_dbConnection = None
_dbConfig = None


def connect(config):
    ''' creates a connection to the database and returns it
    
    Arguments:
    config -- configuration dict

    '''
    global _replace_placeholders
    global _dbtype
    global _dbConnection
    global _dbConfig
    global registrationMetadataTableName
    global origDataTableName
    global hashTableName
    global molblocksTableName
    global conformersTableName
    global lwregSchema

    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    cn = config.get('connection', None)
    if not cn and _dbConnection is not None and _dbConfig == config:
        cn = _dbConnection
    dbtype = _lookupWithDefault(config, 'dbtype').lower()
    if not cn:
        dbnm = config['dbname']
        if dbtype == 'sqlite3':
            uri = False
            if dbnm.startswith('file::'):
                uri = True
            cn = sqlite3.connect(dbnm, uri=uri)
        elif dbtype == "postgresql":
            if psycopg2 is None:
                raise ValueError("psycopg2 package not installed")
            cn = psycopg2.connect(database=dbnm,
                                  host=config.get('host', None),
                                  user=config.get('user', None),
                                  password=config.get('password', None))

    schemaBase = ''
    lwregSchema = ''
    if dbtype == 'postgresql':
        lwregSchema = config.get('lwregSchema', '')
        if lwregSchema:
            schemaBase = config['lwregSchema'] + '.'
    registrationMetadataTableName = schemaBase + _baseregistrationMetadataTableName
    origDataTableName = schemaBase + _baseorigDataTableName
    hashTableName = schemaBase + _basehashTableName
    molblocksTableName = schemaBase + _basemolblocksTableName
    conformersTableName = schemaBase + _baseconformersTableName

    _dbtype = dbtype
    if dbtype == 'postgresql':
        _replace_placeholders = _replace_placeholders_pcts
    else:
        _replace_placeholders = _replace_placeholders_noop
    _dbConnection = cn
    _dbConfig = config
    return cn


def _clear_cached_connection():
    global _dbConnection
    global _dbConfig
    _dbConnection = None
    _dbConfig = None


MolTuple = namedtuple('MolTuple', ('mol', 'datatype', 'rawdata'))


def _process_molblock(molblock, config):
    mol = Chem.MolFromMolBlock(molblock,
                               sanitize=False,
                               removeHs=_lookupWithDefault(config, 'removeHs'))
    return mol


def _parse_mol(mol=None, molfile=None, molblock=None, smiles=None, config={}):
    if mol is not None:
        datatype = 'pkl'
        raw = base64.encodebytes(
            mol.ToBinary(propertyFlags=Chem.PropertyPickleOptions.AllProps))
        raw = raw.decode()
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
    elif isinstance(config, str):
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
        if isinstance(sopt, str):
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
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

    sopts = _get_standardization_list(config)
    for sopt in sopts:
        if isinstance(sopt, str):
            if sopt in standardizationOptions:
                sopt = standardizationOptions[sopt]
            else:
                sopt = getattr(standardization_lib, sopt)()
        mol = sopt(mol)
        if mol is None:
            return None
    return mol


def _get_conformer_hash(mol, numDigits, confId=-1):
    """ returns a hash for a molecule's conformer based on the atomic positions
    
    """
    if not mol.GetNumConformers():
        return ''
    ps = mol.GetConformer(id=confId).GetPositions()
    hps = []
    for p in ps:
        coords = [str(round(x, numDigits)) for x in p]
        hps.append(','.join(coords))
    return sha256(';'.join(sorted(hps)).encode()).hexdigest()


def hash_mol(mol, escape=None, config=None):
    """ returns the hash layers and corresponding hash for the molecule 
    
    Keyword arguments:
    escape -- the escape layer
    config -- configuration dict
    """
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

    layers = RegistrationHash.GetMolLayers(
        mol,
        escape=escape,
        enable_tautomer_hash_v2=_lookupWithDefault(config,
                                                   'useTautomerHashv2'))

    mhash = RegistrationHash.GetMolHash(layers)

    return mhash, layers


def _register_one_conformer(mrn,
                            sMol,
                            molb,
                            cn,
                            curs,
                            config,
                            fail_on_duplicate,
                            confId=-1):
    try:
        chash = _get_conformer_hash(sMol,
                                    _lookupWithDefault(config,
                                                       "numConformerDigits"),
                                    confId=confId)
        regtuple = (mrn, chash, molb)
        qs = '?,?,?'
        if _dbtype != 'postgresql':
            curs.execute(
                _replace_placeholders(
                    f'insert into {conformersTableName} values (NULL,{qs})'),
                regtuple)
            curs.execute('select last_insert_rowid()')
            conf_id = curs.fetchone()[0]
        else:
            # Note that on postgresql the serial ids are increasing, but not sequential
            #  i.e. failed inserts will increment the counter
            curs.execute(
                _replace_placeholders(
                    f'insert into {conformersTableName} values (default,{qs}) returning conf_id'
                ), regtuple)
            conf_id = curs.fetchone()[0]
        cn.commit()
    except _violations:
        cn.rollback()
        if fail_on_duplicate:
            raise
        else:
            curs.execute(
                _replace_placeholders(
                    f'select conf_id from {conformersTableName} where conformer_hash=?'
                ), (chash, ))
            conf_id = curs.fetchone()[0]
    except:
        cn.rollback()
        raise
    return conf_id


def _register_mol(tpl,
                  escape,
                  cn,
                  curs,
                  config,
                  fail_on_duplicate,
                  def_rdkit_version_label=None,
                  def_std_label=None,
                  confId=-1,
                  molCache=None):
    """ does the work of registering one molecule """
    registerConformers = _lookupWithDefault(config, "registerConformers")

    if registerConformers and tpl.mol is not None \
        and not tpl.mol.GetNumConformers():
        raise ValueError(
            "attempt to register a molecule without conformers when registerConformers is set"
        )

    if hasattr(cn, 'autocommit') and cn.autocommit is True:
        logging.warn("setting autocommit on the database connection to False")
        cn.autocommit = False

    standardization_label = _get_standardization_label(config)
    if def_std_label is None:
        curs.execute(
            f"select value from {registrationMetadataTableName} where key='standardization'"
        )
        def_std_label = curs.fetchone()[0]
    if standardization_label == def_std_label:
        standardization_label = None

    if def_rdkit_version_label is None:
        curs.execute(
            f"select value from {registrationMetadataTableName} where key='rdkitVersion'"
        )
        def_rdkit_version_label = curs.fetchone()[0]
    rdkit_version_label = rdkit.__version__
    if rdkit_version_label == def_rdkit_version_label:
        rdkit_version_label = None
    mrn = None
    conf_id = None
    try:
        sMol = standardize_mol(tpl.mol, config=config)
        if confId != -1:
            Chem.AssignStereochemistryFrom3D(sMol, confId)
        if sMol is None:
            return None, None
        try:
            molb = Chem.MolToV3KMolBlock(sMol, confId=confId)
        except KekulizeException:
            return None, None

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
                    f'insert into {hashTableName} values (NULL,{qs})'),
                regtuple)
            curs.execute('select last_insert_rowid()')
            mrn = curs.fetchone()[0]
        else:
            # Note that on postgresql the serial ids are increasing, but not sequential
            #  i.e. failed inserts will increment the counter
            curs.execute(
                _replace_placeholders(
                    f'insert into {hashTableName} values (default,{qs}) returning molregno'
                ), regtuple)
            mrn = curs.fetchone()[0]

        curs.execute(
            _replace_placeholders(
                f'insert into {origDataTableName} (molregno, data, datatype) values (?, ?, ?)'
            ), (mrn, tpl.rawdata, tpl.datatype))
        curs.execute(
            _replace_placeholders(
                f'insert into {molblocksTableName} values (?, ?, ?)'),
            (mrn, molb, standardization_label))

        cn.commit()
    except _violations:
        cn.rollback()
        if fail_on_duplicate and not (registerConformers
                                      and sMol.GetNumConformers()):
            raise
        else:
            curs.execute(
                _replace_placeholders(
                    f'select molregno from {hashTableName} where fullhash=?'),
                (mhash, ))
            mrn = curs.fetchone()[0]
    except:
        cn.rollback()
        raise

    if molCache is not None:
        molCache.append(sMol)
    if registerConformers and sMol.GetNumConformers():
        conf_id = _register_one_conformer(mrn,
                                          sMol,
                                          molb,
                                          cn,
                                          curs,
                                          config,
                                          fail_on_duplicate,
                                          confId=confId)

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
             confId=-1,
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
    fail_on_duplicate -- if true then an exception is raised when trying to register a duplicate
    confId     -- the conformer ID to use when in registerConformers mode
    no_verbose -- if this is False then the registry number will be printed
    """
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

    tpl = _parse_mol(mol=mol,
                     molfile=molfile,
                     molblock=molblock,
                     smiles=smiles,
                     config=config)

    cn = connect(config)
    curs = cn.cursor()
    if tpl.mol is None:
        return RegistrationFailureReasons.PARSE_FAILURE

    mrn, conf_id = _register_mol(tpl,
                                 escape,
                                 cn,
                                 curs,
                                 config,
                                 fail_on_duplicate,
                                 confId=confId)
    if mrn is None:
        return RegistrationFailureReasons.FILTERED
    if not _lookupWithDefault(config, "registerConformers"):
        res = mrn
    else:
        res = mrn, conf_id
    if not no_verbose:
        print(res)
    return res


def register_multiple_conformers(config=None,
                                 mol=None,
                                 escape=None,
                                 fail_on_duplicate=True,
                                 no_verbose=True):
    """ registers all of the conformers of a multi-conformer molecule
    Using this function only makes sense when registerConformers is enabled.
    
    Keyword arguments:
    config     -- configuration dict
    mol        -- RDKit molecule object (must have at least one conformer)
    escape     -- the escape layer
    fail_on_duplicate -- if true then an exception is raised when trying to register a duplicate
    no_verbose -- if this is False then the registry number will be printed
    """
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)
    if not _lookupWithDefault(config, "registerConformers"):
        raise ValueError(
            'register_multiple_conformers can only be used when registerConformers is enabled'
        )

    tpl = _parse_mol(mol=mol, config=config)
    if tpl.mol is None:
        return RegistrationFailureReasons.PARSE_FAILURE

    cn = connect(config)
    curs = cn.cursor()

    # start by registering the first conformer in order to
    # get the molregno that we'll use later
    rc = []
    mrns = {}
    confsDone = set()
    confMrns = []
    res = []
    for i, conf in enumerate(tpl.mol.GetConformers()):
        Chem.AssignStereochemistryFrom3D(tpl.mol, conf.GetId())
        smi = Chem.MolToSmiles(tpl.mol)
        if smi not in mrns:
            mrn, conf_id = _register_mol(tpl,
                                         escape,
                                         cn,
                                         curs,
                                         config,
                                         fail_on_duplicate,
                                         confId=conf.GetId(),
                                         molCache=rc)
            if mrn is not None:
                mrns[smi] = mrn
                confsDone.add(i)
                res.append((mrn, conf_id))
        else:
            mrn = mrns[smi]
        confMrns.append(mrn)
    if not len(res):
        return RegistrationFailureReasons.FILTERED

    sMol = rc[0]
    for i, conf in enumerate(sMol.GetConformers()):
        if i in confsDone:
            # we already registered the first conformer
            continue
        molb = Chem.MolToV3KMolBlock(sMol, confId=conf.GetId())
        mrn = confMrns[i]
        conf_id = _register_one_conformer(mrn,
                                          sMol,
                                          molb,
                                          cn,
                                          curs,
                                          config,
                                          fail_on_duplicate,
                                          confId=conf.GetId())
        res.append((mrn, conf_id))

    if not no_verbose:
        print(res)
    return tuple(res)


def bulk_register(config=None,
                  mols=None,
                  sdfile=None,
                  smilesfile=None,
                  escape_property=None,
                  fail_on_duplicate=True,
                  no_verbose=True,
                  show_progress=False):
    """ registers multiple new molecules, assuming they don't already exist,
    and returns the new registry numbers (molregno)
    
    The result tuple includes a single entry for each molecule in the input.
    That entry can be one of the following:
      - the registry number (molregno) of the registered molecule
      - RegistrationFailureReasons.DUPLICATE if fail_on_duplicate is True and a
        molecule is a duplicate
      - RegistrationFalureReasons.PARSE_FAILURE if there was a problem processing
        the molecule 
    
    only one of the molecule format objects should be provided

    Keyword arguments:
    config         -- configuration dict
    mols           -- an iterable of RDKit molecule objects
    sdfile         -- SDF filename
    escape_property -- the molecule property to use as the escape layer
    fail_on_duplicate -- if true then RegistraionFailureReasons.DUPLICATE will be returned 
                       for each already-registered molecule, otherwise the already existing
                       structure ID will be returned
    no_verbose     -- if this is False then the registry numbers will be printed
    show_progress   -- if this is True then a progress bar will be shown for the molecules
    """

    if mols:
        pass
    elif sdfile:
        mols = Chem.ForwardSDMolSupplier(sdfile,
                                         removeHs=_lookupWithDefault(
                                             config, 'removeHs'))
    elif smilesfile:
        mols = _get_mols_from_smilesfile(smilesfile)
    else:
        raise ValueError('No input molecules provided!')

    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

    res = []
    cn = connect(config)
    curs = cn.cursor()

    curs.execute(
        f"select value from {registrationMetadataTableName} where key='standardization'"
    )
    def_std_label = curs.fetchone()[0]
    curs.execute(
        f"select value from {registrationMetadataTableName} where key='rdkitVersion'"
    )
    def_rdkit_version_label = curs.fetchone()[0]
    for mol in tqdm(mols, disable=not show_progress):
        if mol is None:
            res.append(RegistrationFailureReasons.PARSE_FAILURE)
            continue
        tpl = _parse_mol(mol=mol, config=config)
        try:
            if escape_property is not None and mol.HasProp(escape_property):
                escape = mol.GetProp(escape_property)
            else:
                escape = None
            mrn, conf_id = _register_mol(
                tpl,
                escape,
                cn,
                curs,
                config,
                fail_on_duplicate,
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

    cn = connect(config)
    curs = cn.cursor()

    qs = _replace_placeholders(','.join('?' * len(ids)))
    curs.execute(
        f'select molregno,conf_id from {conformersTableName} where molregno in ({qs})',
        ids)
    res = curs.fetchall()
    return res


def registration_counts(config=None):
    """ returns the number of entries in the registration database

    the result is the number of molecules if registerConformers is not set,
    and a tuple of (number of molecules, number of conformers) if it is

    Keyword arguments:
    config     -- configuration dict
    """
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    cn = connect(config)
    curs = cn.cursor()
    curs.execute(f'select count(*) from {hashTableName}')
    nHashes = curs.fetchone()[0]

    if _lookupWithDefault(config, "registerConformers"):
        curs.execute(f'select count(*) from {conformersTableName}')
        nConfs = curs.fetchone()[0]
        return nHashes, nConfs
    else:
        return nHashes


def get_all_registry_numbers(config=None):
    """ returns a tuple with all of the registry numbers (molregnos) in the database
        
    Keyword arguments:
    config     -- configuration dict
    """
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    cn = connect(config)
    curs = cn.cursor()
    curs.execute(f'select molregno from {hashTableName}')
    res = curs.fetchall()
    return tuple(sorted(x[0] for x in res))


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
    and returns the corresponding registry numbers
        
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
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

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

        cn = connect(config)
        curs = cn.cursor()
        if not _lookupWithDefault(config, "registerConformers") or \
           not sMol.GetNumConformers():
            mhash, hlayers = hash_mol(sMol, escape=escape, config=config)
            if layers == 'ALL':
                queryText = _replace_placeholders(
                    f'select molregno from {hashTableName} where fullhash=?')
                queryVals = (mhash, )
            else:
                vals = []
                query = []
                if isinstance(layers, str):
                    layers = layers.upper().split(',')
                if not hasattr(layers, '__len__'):
                    layers = [layers]
                for lyr in layers:
                    if isinstance(lyr, str):
                        k = getattr(HashLayer, lyr)
                    else:
                        k = lyr
                        lyr = str(lyr).split('.')[-1]
                    vals.append(hlayers[k])
                    query.append(f'"{lyr}"=?')

                query = _replace_placeholders(' and '.join(query))
                queryText = f'select molregno from {hashTableName} where {query}'
                queryVals = vals
            curs.execute(queryText, queryVals)
            res = [x[0] for x in curs.fetchall()]
        else:
            # do a conformer query
            chash = _get_conformer_hash(
                sMol, _lookupWithDefault(config, "numConformerDigits"))
            curs.execute(
                _replace_placeholders(
                    f'select molregno,conf_id from {conformersTableName} where conformer_hash=?'
                ), (chash, ))
            res = [tuple(x) for x in curs.fetchall()]

    if not no_verbose:
        if res:
            print(' '.join(str(x) for x in res))
        else:
            print('not found')

    return res


def _parsePickleFromDB(data):
    # there was a bug in early versions of the code that stored the binary data really badly
    if data[:2] == r'\x':
        byts = []
        for i in range(2, len(data), 2):
            byts.append(int(data[i:i + 2], base=16))
        data = bytes(byts)
    else:
        data = base64.decodebytes(data.encode())
    return data


def retrieve(config=None,
             ids=None,
             id=None,
             as_submitted=False,
             as_hashes=False,
             no_verbose=True):
    """ returns the molecule data for one or more registry ids (molregnos)
    The return value is a dictionary of (data, format) 2-tuples with molregnos as keys

    only one of id or ids should be provided

    If registerConformers is set the conformers can be retrieved by providing
    the tuples of (molregno, conf_id) and the return value will be a dictionary
    of (data, 'mol') 2-tuples with (molregno, conf_id) tuples as keys

    Keyword arguments:
    config       -- configuration dict
    ids          -- an iterable of registry ids (molregnos)
    id           -- a registry id (molregno)
    as_submitted -- if True, then the structure will be returned as registered
    as_hashes    -- if True, then the hashes will be returned (as a dict) instead of the structures
    no_verbose   -- if this is False then the registry number will be printed
    """
    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

    registerConformers = _lookupWithDefault(config, "registerConformers")

    if id is not None:
        if registerConformers:
            try:
                ids = [(int(id[0]), int(id[1]))]
                getConfs = True
            except TypeError:
                ids = [int(id)]
                getConfs = False
        else:
            ids = [int(id)]
            getConfs = False
    elif ids is not None:
        if registerConformers:
            try:
                ids = [(int(x), int(y)) for x, y in ids]
                getConfs = True
            except TypeError:
                ids = [int(x) for x in ids]
                getConfs = False
        else:
            if isinstance(ids, str):
                ids = [int(x) for x in ids.split(',')]
            getConfs = False

    cn = connect(config)
    curs = cn.cursor()
    if not getConfs:
        if as_hashes:
            qry = f'* from {hashTableName}'
        elif as_submitted:
            qry = f'molregno,data,datatype from {origDataTableName}'
        else:
            qry = f"molregno,molblock,'mol' from {molblocksTableName}"
        qs = _replace_placeholders(','.join('?' * len(ids)))
        curs.execute(f'select {qry} where molregno in ({qs})', ids)
    else:
        qry = f"molregno,conf_id,molblock from {conformersTableName}"
        qs = _replace_placeholders(','.join('?' * len(ids)))
        curs.execute(
            f'select {qry} where molregno in ({qs}) and conf_id in ({qs})',
            ([x for x, y in ids] + [y for x, y in ids]))

    res = curs.fetchall()
    if not no_verbose:
        if res:
            for entry in res:
                print(entry)
        else:
            print('not found')
    resDict = {}
    if as_hashes:
        colns = [x[0] for x in curs.description]
        for row in res:
            rowd = {}
            for i, coln in enumerate(colns):
                if coln == 'rdkitversion':
                    continue
                if coln == 'molregno':
                    mrn = row[i]
                    continue
                if row[i] is None:
                    continue
                rowd[coln] = row[i]
            resDict[mrn] = rowd
    else:
        if not getConfs:
            for mrn, data, fmt in res:
                if as_submitted and fmt == 'pkl':
                    data = _parsePickleFromDB(data)
                resDict[mrn] = (data, fmt)
        else:
            for mrn, confId, molb in res:
                resDict[(mrn, confId)] = (molb, 'mol')

    return resDict


def _registerMetadata(curs, config):
    dc = defaultConfig()
    dc.update(config)
    dc['standardization'] = _get_standardization_label(dc)
    for k, v in dc.items():
        if k in ('connection', 'user', 'password'):
            continue
        curs.execute(
            _replace_placeholders(
                f'insert into {registrationMetadataTableName} values (?,?)'),
            (str(k), str(v)))

    curs.execute(
        _replace_placeholders(
            f'insert into {registrationMetadataTableName} values (?,?)'),
        ('rdkitVersion', rdkit.__version__))


def _initdb(config=None, confirm=False):
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
    elif isinstance(config, str):
        config = _configure(filename=config)

    _check_config(config)

    cn = connect(config)
    curs = cn.cursor()

    if lwregSchema and _dbtype == 'postgresql':
        curs.execute(f'create schema if not exists {lwregSchema}')
    curs.execute(f'drop table if exists {registrationMetadataTableName}')
    curs.execute(
        f'create table {registrationMetadataTableName} (key text, value text)')
    _registerMetadata(curs, config)
    cn.commit()

    if _dbtype != 'postgresql':
        curs.execute(f'drop table if exists {hashTableName}')
    else:
        curs.execute(f'drop table if exists {hashTableName} cascade')
    if _dbtype != 'postgresql':
        curs.execute(
            f'''create table {hashTableName} (molregno integer primary key, fullhash text unique, 
            formula text, canonical_smiles text, no_stereo_smiles text, 
            tautomer_hash text, no_stereo_tautomer_hash text, "escape" text, sgroup_data text, rdkitVersion text)'''
        )
        curs.execute(
            f'''create unique index {hashTableName}_fullhash_idx on {hashTableName} 
                (fullhash)''')
    else:
        curs.execute(
            f'''create table {hashTableName} (molregno serial primary key, fullhash text unique, 
            formula text, canonical_smiles text, no_stereo_smiles text, 
            tautomer_hash text, no_stereo_tautomer_hash text, "escape" text, sgroup_data text, rdkitVersion text)'''
        )
        curs.execute(
            f'''create index {hashTableName.replace(".","_")}_fullhash_idx on {hashTableName} 
                using hash(fullhash)''')
    curs.execute(f'drop table if exists {origDataTableName}')
    if _dbtype != 'postgresql':
        curs.execute(
            f'create table {origDataTableName} (molregno integer unique not null, data text, datatype text, timestamp DATETIME DEFAULT CURRENT_TIMESTAMP, foreign key(molregno) references {hashTableName} (molregno))'
        )
    else:
        curs.execute(
            f'create table {origDataTableName} (molregno integer unique not null references {hashTableName} (molregno), data text, datatype text, timestamp TIMESTAMP default now())'
        )
    curs.execute(f'drop table if exists {molblocksTableName}')
    curs.execute(
        f'create table {molblocksTableName} (molregno integer unique not null references {hashTableName} (molregno), molblock text, standardization text, foreign key(molregno) references {hashTableName} (molregno))'
    )

    curs.execute(f'drop table if exists {conformersTableName}')
    if _lookupWithDefault(config, "registerConformers"):
        if _dbtype != 'postgresql':
            curs.execute(
                f'''create table {conformersTableName} (conf_id integer primary key, molregno integer not null, 
                   conformer_hash text not null unique, molblock text, foreign key(molregno) references {hashTableName} (molregno))'''
            )
            curs.execute(
                f'''create unique index {conformersTableName}_hash_idx on {conformersTableName} 
                    (conformer_hash)''')
        else:
            curs.execute(
                f'''create table {conformersTableName} (conf_id serial primary key, molregno integer references {hashTableName} (molregno), 
                   conformer_hash text not null unique, molblock text)''')
            curs.execute(
                f'''create index {conformersTableName.replace(".","_")}_hash_idx on {conformersTableName} 
                    using hash(conformer_hash)''')

    cn.commit()
    return True


def initdb(config=None):
    """ initializes the registration database    

    NOTE you will be prompted to confirm this action since this call can destroy any 
    existing information in the registration database

    Keyword arguments:
    config  -- configuration dict
    """
    print(
        "This will destroy any existing information in the registration database."
    )
    response = input("  are you sure? [yes/no]: ")
    if response == 'yes':
        return _initdb(config=config, confirm=True)
    else:
        print("cancelled")
        return False


def _check_config(config):
    ''' checks that the configuration is valid and no forbidden combinations are present
        
        Selected options in the configuration dict are checked against a list of
        forbidden combinations. If any of the combinations are present a ValueError
        is raised.

    '''

    if not config:
        config = _configure()
    elif isinstance(config, str):
        config = _configure(filename=config)

    if config.get("dbtype", "sqlite3") not in ('sqlite3', 'postgresql'):
        raise ValueError(
            "Possible values for dbtype are sqlite3 and postgresql")
