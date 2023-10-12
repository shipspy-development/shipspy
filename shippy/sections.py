import numpy as np
import xarray as xr
import pandas as pd
import argparse
#from .helpers import get_args

def configure_sections_parser(parser):
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Postprocessed file (netCDF)", default=None,
                        required=True)

    parser.add_argument("-o", "--outputfile", metavar="OUTPUT_FILE",
                        help="Output file (netCDF)", default=None,
                        required=True)
    
    parser.add_argument("-s", "--sectionfile", metavar="SECTION_FILE",
                        help="txt file with campaign and section information", default=None,
                        required=True)
    
    parser.add_argument("-d", "--dimensioname", metavar="start_time",
                        help="Name of dimension along which section is assigned", default='time',
                        required=False)

    parser.add_argument('-v', '--verbose', metavar="DEBUG",
                        help='Set the level of verbosity [DEBUG, INFO,'
                        ' WARNING, ERROR]',
                        required=False, default="INFO")

    parser.set_defaults(func=_run)


def _run(args):
    infilename = args.inputfile
    outfilename = args.outputfile
    secfilename = args.sectionfile
    dim_name = args.dimensioname
    print(dim_name)