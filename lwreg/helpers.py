from . import standardization_lib
from . import utils
import yaml
import getpass
import json

def interactive_config():
    config = {}
    config['dbname'] = input("Enter the name of your database: \n")
    dbtype_option = int(input("Choose your database type: \n options: \n 1: sqlite \n 2: postgresql \n"))
    if dbtype_option == 1:
        config["dbtype"] = "sqlite3"
    elif dbtype_option == 2:
        config["dbtype"] = "postgresql"
    else:
        raise ValueError('Choosen option is invalid')
    std_option = int(input("Choose your standardization: \n options: 1: none, 2: sanitize, 3: fragment, 4: charge, 5: tautomer, 6: super \n"))
    std_option_mapping = {1: "none",
                          2: "sanitize",
                          3: "fragment",
                          4: "charge",
                          5: "tautomer",
                          6: "super"}
    if std_option >= 1 and std_option <= 6:
        config["standardization"] = std_option_mapping[std_option]
    else:
        raise ValueError('Choosen option is invalid')
    Hs_option = int(input("Choose if you want to remove Hs: \n options: \n 1: no, 2: yes \n"))
    if Hs_option == 1:
        config["removeHs"] = 0
    elif Hs_option == 2:
        config["removeHs"] = 1
    else:
        raise ValueError('Choosen option is invalid')
    tautomerhash_option = int(input("Do you want to use the TautomerHashv2? \n options \n 1: no, 2: yes \n"))
    if tautomerhash_option == 1:
        config["useTautomerHashv2"] = 0
    elif tautomerhash_option == 2:
        config["useTautomerHashv2"] = 1
    else:
        raise ValueError('Choosen option is invalid')
    conf_option = int(input("Do you want to register conformers? \n options: 1: no, 2: yes \n"))
    if conf_option == 1:
        config["registerConformers"] = 0
    elif conf_option == 2:
        config["registerConformers"] = 1
    else: 
        raise ValueError('Choosen option is invalid')
    numConfDigits_opt = int(input("Do you want to set the number of conformer digits? \n options: \n 1: no (keep default 3), 2: yes \n"))
    if numConfDigits_opt == 1:
        config["numConformerDigits"] = 3
    elif numConfDigits_opt == 2:
        config["numConformerDigits"] = int(input("Please enter the number of digits you want to keep \n"))
    else:
        raise ValueError('Choosen option is invalid')
    if config["registerConformers"] != 1:
        hashConf_opt = int(input("Do you want the conformer to be part of the basic identity hash? \n options: \n 1: no, 2: yes \n"))
        if hashConf_opt == 1:
            config["hashConformer"] = 0
        elif hashConf_opt == 2:
            config["hashConformer"] = 1
            raise ValueError('Choosen option is invalid')
    else:
        config["hashConformer"] = 0
    if config["dbtype"] == "postgresql":
        config["lwregSchema"] = "public"
        newName = input("If you would like lwreg to use its own schema (default: public), enter it here: \n")
        if newName:
            config["lwregSchema"] = newName
    config["host"] = input("Please enter the name of your database host machine \n")
    config["user"] = input("Please enter your username \n")
    print("Please enter your password \n")
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