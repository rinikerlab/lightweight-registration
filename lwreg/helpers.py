from . import utils
import json

def interactive_config():
    """ an interactive configuration assistant returning a valid configuration dictionary for an lwreg instance
    """
    print("This is an interactive configuration assisstant for lwreg. \n Please manually add host, user and password to your configuration.")
    config = utils._defaultConfig
    config['dbname'] = input("Enter the name of your database: ")
    dbtype_option = input("Choose your database type: ([sqlite3], postgresql) ")
    if not dbtype_option:
        config["dbtype"] = "sqlite3"
    elif dbtype_option in ["sqlite3","postgresql"]:
        config["dbtype"] = dbtype_option
    else:
        raise ValueError('Selected option is invalid')
    std_option = input("Choose your standardization: (none, sanitize, [fragment], charge, tautomer, super) ")
    if not std_option:
        config["standardization"] = "fragment"
    elif std_option in ["none", "sanitize", "fragment", "charge", "tautomer", "super"]:
        config["standardization"] = std_option
    else:
        raise ValueError('Selected option is invalid')
    Hs_option = input("Do you want to remove Hs? ([yes]/no) ")
    if not Hs_option:
        config["removeHs"] = 1
    elif Hs_option == "no":
        config["removeHs"] = 0
    elif Hs_option == "yes":
        config["removeHs"] = 1
    else:
        raise ValueError('Selected option is invalid')
    tautomerhash_option = input("Do you want to use the TautomerHashv2? (yes/[no]) ")
    if not tautomerhash_option:
        config["useTautomerHashv2"] = 0
    elif tautomerhash_option == "no":
        config["useTautomerHashv2"] = 0
    elif tautomerhash_option == "yes":
        config["useTautomerHashv2"] = 1
    else:
        raise ValueError('Selected option is invalid')
    conf_option = input("Do you want to register conformers? (yes/[no]) ")
    if not conf_option:
        config["registerConformers"] = 0
    elif conf_option == "no":
        config["registerConformers"] = 0
    elif conf_option == "yes":
        config["registerConformers"] = 1
        numConfDigits_opt = input("How many conformer digits would you like to include? ([3]) ")
        if not numConfDigits_opt:
            config["numConformerDigits"] = 3
        else:
            config["numConformerDigits"] = int(numConfDigits_opt)
    else: 
        raise ValueError('Selected option is invalid')

    if config["dbtype"] == "postgresql":
        config["lwregSchema"] = ""
        newName = input("If you would like lwreg to use its own schema, enter it here: ")
        if newName:
            config["lwregSchema"] = newName
    utils._check_config(config)
    
    return config

def write_configfile(config,config_filename="config.json"):
    """ Writes a stripped version of a configuration dictionary to a json file.
    The information in the file will allow the user to connect to the lwreg instance after addding username and password. 
    From the lwreg instance itself, all initial configuration information can be retrieved.
    """
    config_stripped = {}
    config_stripped["dbname"] = config["dbname"]
    config_stripped["dbtype"] = config["dbtype"]
    config_stripped["host"] = config["host"]
    config_stripped["lwregSchema"] = config["lwregSchema"]

    with open(config_filename,"w") as f:
        json.dump(config_stripped,f)

def load_configfile(config_filename):
    """ Loads a configuration dictionary from a json file.
    """
    with open(config_filename,"r") as f:
        config =  json.load(f)
    return config