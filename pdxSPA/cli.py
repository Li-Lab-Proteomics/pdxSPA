import os
import sys
import re
import click

from . import _version
from .pdxSPA import SPA, X_SPA, S_SPA

__author__ = 'Xi Cheng'
__copyright__ = _version.__copyright__
__credits__ = ['Xi Cheng']
__license__ = 'SJTU'
__version__ = _version.__version__
__maintainer__ = 'Xi Cheng'
__email__ = 'chengxi0237@sjtu.edu.cn'

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(name='main')
@click.version_option(__version__)
def cli():
    pass

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--input_file', '-i', type=click.Path(exists=True),
              default=None, show_default=True, 
              help='Pass in the evidence file after quantifile PDX samples by MaxQuant.')
@click.option('--outdir', '-o', type=click.Path(exists=True),
              default='.', show_default=True, 
              help='Pass in the output directory.')

@click.option('--method', '-m', type=int,
              default=0, show_default=True, 
              help='Choose methods between 0 (SPA), 1 (X-SPA) and 2 (S-SPA).')
@click.version_option(__version__)

def assignSharedPeptides(input_file, outdir, method):
    if input_file:
        try:
            f = open(input_file)
            f.close()
        except FileNotFoundError:
            print("Input file is not found.")
        except PersmissionError:
            print("Do not have permission to access the input file.")
    if outdir[-1] != '/':
        outdir += '/'

    if method == 0:
        SPA(input_file, outdir)
    elif method == 1:
        X_SPA(input_file, outdir)
    elif method == 2:
        S_SPA(input_file, outdir)
    else:
        raise ValueError('Invalid method')

def main():
    assignSharedPeptides()
