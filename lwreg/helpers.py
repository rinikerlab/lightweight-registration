from . import standardization_lib
from . import utils
import yaml
import getpass
import json

def interactive_config():
    config = {}
    print("Enter the name of your database: ")
    config['dbname'] = input()
    print("Choose your database type: \n options: \n 1: sqlite \n 2: postgresql")
    dbtype_option = int(input())
    if dbtype_option == 1:
        config["dbtype"] = "sqlite3"
    elif dbtype_option == 2:
        config["dbtype"] = "postgresql"
    else:
        raise ValueError('Choosen option is invalid')
    print("Choose your standardization: \n options: 1: none, 2: sanitize, 3: fragment, 4: charge, 5: tautomer, 6: super")
    std_option = int(input())
    if std_option == 1:
        config["standardization"] = 'none' # standardization_lib.NoStandardization()
    elif std_option == 2:
        config["standardization"] = 'sanitize' # standardization_lib.RDKitSanitize()
    elif std_option == 3:
        config["standardization"] = 'fragment' # standardization_lib.FragmentParent()
    elif std_option == 4:
        config["standardization"] = 'charge' # standardization_lib.ChargeParent()
    elif std_option == 5:
        config["standardization"] = 'tautomer' # standardization_lib.TautomerParent()
    elif std_option == 6:
        config["standardization"] = 'super' # standardization_lib.SuperParent()
    else:
        raise ValueError('Choosen option is invalid')
    print("Choose if you want to remove Hs: \n options: \n 1: no, 2: yes")
    Hs_option = int(input())
    if Hs_option == 1:
        config["removeHs"] = 0
    elif Hs_option == 2:
        config["removeHs"] = 1
    else:
        raise ValueError('Choosen option is invalid')
    print("Do you want to use the TautomerHashv2? \n options \n 1: no, 2: yes")
    tautomerhash_option = int(input())
    if tautomerhash_option == 1:
        config["useTautomerHashv2"] = 0
    elif tautomerhash_option == 2:
        config["useTautomerHashv2"] = 1
    else:
        raise ValueError('Choosen option is invalid')
    print("Do you want to register conformers? \n options: 1: no, 2: yes")
    conf_option = int(input())
    if conf_option == 1:
        config["registerConformers"] = 0
    elif conf_option == 2:
        config["registerConformers"] = 1
    else: 
        raise ValueError('Choosen option is invalid')
    print("Do you want to set the number of conformer digits? \n options: \n 1: no (keep default 3), 2: yes ")
    numConfDigits_opt = int(input())
    if numConfDigits_opt == 1:
        config["numConformerDigits"] = 3
    elif numConfDigits_opt == 2:
        print("Please enter the number of digits you want to keep")
        config["numConformerDigits"] = int(input())
    else:
        raise ValueError('Choosen option is invalid')
    if config["registerConformers"] != 1:
        print("Do you want the conformer to be part of the basic identity hash? \n options: \n 1: no, 2: yes")
        hashConf_opt = int(input())
        if hashConf_opt == 1:
            config["hashConformer"] = 0
        elif hashConf_opt == 2:
            config["hashConformer"] = 1
            raise ValueError('Choosen option is invalid')
    else:
        config["hashConformer"] = 0
    if config["dbtype"] == "postgresql":
        print("Do you want to place lwreg in its own schema? \n options: \n 1: no (default: public), 2: yes")
        schema_opt = int(input())
        if schema_opt == 1:
            config["lwregSchema"] = "public"
        elif schema_opt == 2:
            print("Please enter your chosen schema name")
            config["lwregSchema"] = input()
        else:
            raise ValueError('Choosen option is invalid')
    print("Please enter the name of your database host machine")
    config["host"] = input()
    print("Please enter your username")
    config["user"] = input()
    print("Please enter your password")
    config["password"] = getpass.getpass()
    utils._check_config(config)
    return config

def write_configfile(config,config_filename="lwreg_config.yml"):
    config_stripped = {}
    config_stripped["dbname"] = config["dbname"]
    config_stripped["dbtype"] = config["dbtype"]
    config_stripped["host"] = config["host"]
    config_stripped["user"] = config["user"]
    config_stripped["lwregSchema"] = config["lwregSchema"]

    with open(config_filename,"w") as f:
        yaml.dump(config_stripped,f,default_flow_style=False)
    return

def load_configfile(config_filename):
    with open(config_filename,"r") as file:
        config =  yaml.safe_load(file)
    return config