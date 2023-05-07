#! /bin/env python

# Copyright (C) 2022 Greg Landrum
# All rights reserved
# This file is part of lwreg.
# The contents are covered by the terms of the MIT license
# which is included in the file LICENSE,
import click
import logging
import sys
from . import utils


@click.group()
@click.option(
    "--config",
    default='',
)
def cli(config=''):
    utils._configure(filename=config)


@cli.command()
@click.option(
    "--confirm",
    default='no',
)
def initdb(confirm='no'):
    if confirm != 'yes':
        click.echo("initdb not confirmed, aborting")
        return
    utils.initdb(confirm=True)


@cli.command()
@click.option(
    "--ids",
    default=None,
)
@click.option(
    "--id",
    default=None,
)
@click.option("--as_submitted", default=False, is_flag=True)
@click.option("--no-verbose", default=False, is_flag=True)
def retrieve(**kwargs):
    return utils.retrieve(**kwargs)


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
@click.option("--no-verbose", default=False, is_flag=True)
def query(**kwargs):
    return utils.query(**kwargs)


@cli.command()
@click.option(
    "--smiles",
    default=None,
)
@click.option("--fail-on-duplicate", default=True, is_flag=True)
@click.option(
    "--escape",
    default=None,
)
@click.option("--no-verbose", default=False, is_flag=True)
def register(**kwargs):
    try:
        return utils.register(**kwargs)
    except utils._violations:
        logging.error("Compound already registered")
        sys.exit(1)


@cli.command()
@click.option('--who', default='world')
def greet(who):
    print(f'hello {who}')


if __name__ == '__main__':
    cli()