# Installing lwreg

## Clone this github repo on your computer

We don't have an easy-to-install package yet, so you have to clone this repo locally on your computer. We'll assume you already know how to that since introducing git is beyond the scope of this document. ;-)


## Installing conda and mamba

We highly recommend setting up a conda environment in which to run lwreg. Here's how to do that:

### Install mambaforge

If you already have a conda distribution installed on your machine, you can skip this step.
Otherwise follow the instructions for installing mambaforge here:
https://mamba.readthedocs.io/en/latest/installation.html

### Install mamba (optional, but recommended)

If you don't already have mamba installed (you will if you installed mambaforge), we really recommend doing so; it will allow you to create the environment and install new packages much, much more quickly.

You can do this as follows:
```
conda install -n base -c conda-forge mamba
```
*Note*: the mamba page recommends that you not do this, but I (Greg) have never had a problem with it. Your mileage may vary.

## Create the required environment for lwreg

Open a shell where you can run conda or mamba, go to the directory where you have the clone of this repo, and do:
```
mamba env create -n py311_lwreg --file environment.yml
```
This will take a while, but it will install everything you need to use lwreg, work through the examples in the documentation, and use the demo jupyter notebook.

## Installing lwreg itself.

Again in the shell in the directory with the clone of the repo, activate the environment you just created:
```
conda activate py311_lwreg
```
Then run this command in this directory:
```
pip install --editable .
```

## Verifying that everything worked.

At this point you should be able to type:
```
lwreg --help
```
and see the documentation for the command-line interface to lwreg.

Now you're ready to go... the `README.md` file has more documentation on how to get started with lwreg.
