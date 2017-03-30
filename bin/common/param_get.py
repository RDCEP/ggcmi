#!/usr/bin/env python
# This is used by the psims shell script to extract values from a params file
#

import argparse
import re
import sys
import ruamel.yaml

def get_value(dictionary, key):
    try:
        if isinstance(dictionary[key], list):
            return " ".join(dictionary[key])
        else:
            return dictionary[key]
    except KeyError:
         return ""

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Input params file')
parser.add_argument('--key', help='Key to replace')
options    = parser.parse_args()
inputfile  = options.input
key        = options.key

# Load yaml
inputfile = open(inputfile)
params = ruamel.yaml.load(inputfile, ruamel.yaml.RoundTripLoader)
inputfile.close()

# Handle actions
if ':' in key:
    try:
        (k,v) = key.split(':')
        print get_value(params[k], v)
    except KeyError:
        print ""
else:
    print get_value(params, options.key)
