import numpy as np
import xarray as xr
import pandas as pd
import argparse
#from .helpers import get_args

def configure_dship_parser(parser):
    parser.add_argument("-i", "--inputfile", metavar="INPUT_FILE",
                        help="Postprocessed file (netCDF)", default=None,
                        required=True)

    parser.add_argument("-o", "--outputfile", metavar="OUTPUT_FILE",
                        help="Output file (netCDF)", default=None,
                        required=True)
    

    parser.add_argument('-v', '--verbose', metavar="DEBUG",
                        help='Set the level of verbosity [DEBUG, INFO,'
                        ' WARNING, ERROR]',
                        required=False, default="INFO")

def get_args():
    parser = argparse.ArgumentParser()
    configure_dship_parser(parser)

    parsed_args = vars(parser.parse_args())

    return parsed_args
    
def main():
    args = get_args()
    
    infilename = args['inputfile']
    print(infilename)

if __name__ == "__main__":
    exit(main())
