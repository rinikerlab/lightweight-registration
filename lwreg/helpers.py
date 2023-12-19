from . import standardization_lib
from . import utils
import yaml
import getpass
import json

def interactive_config():
    config = utils._defaultConfig
    config['dbname'] = input("Enter the name of your database: ")
    dbtype_option = int(input("Choose your database type: ([1:sqlite], 2:postgresql) "))
    if dbtype_option == 1:
        config["dbtype"] = "sqlite3"
    elif dbtype_option == 2:
        config["dbtype"] = "postgresql"
    else:
        raise ValueError('Choosen option is invalid')
    std_option = int(input("Choose your standardization: (1:none, 2:sanitize, [3:fragment], 4:charge, 5:tautomer, 6:super) "))
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
    Hs_option = input("Do you want to remove Hs? ([Yes]/no) ")
    if Hs_option == "no":
        config["removeHs"] = 0
    elif Hs_option == "Yes":
        config["removeHs"] = 1
    else:
        raise ValueError('Choosen option is invalid')
    tautomerhash_option = input("Do you want to use the TautomerHashv2? (Yes/[no]) ")
    if tautomerhash_option == "no":
        config["useTautomerHashv2"] = 0
    elif tautomerhash_option == "Yes":
        config["useTautomerHashv2"] = 1
    else:
        raise ValueError('Choosen option is invalid')
    conf_option = input("Do you want to register conformers? (Yes/[no]) ")
    if conf_option == "no":
        config["registerConformers"] = 0
    elif conf_option == "Yes":
        config["registerConformers"] = 1
    else: 
        raise ValueError('Choosen option is invalid')
    numConfDigits_opt = input("Do you want to set the number of conformer digits? (Yes/[no]) ")
    if numConfDigits_opt == "no":
        config["numConformerDigits"] = 3
    elif numConfDigits_opt == "Yes":
        config["numConformerDigits"] = int(input("Please enter the number of digits you want to keep "))
    else:
        raise ValueError('Choosen option is invalid')
    if config["dbtype"] == "postgresql":
        config["lwregSchema"] = ""
        newName = input("If you would like lwreg to use its own schema, enter it here: ")
        if newName:
            config["lwregSchema"] = newName
    config["host"] = input("Please enter the name of your database host machine ")
    config["user"] = input("Please enter your username ")
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