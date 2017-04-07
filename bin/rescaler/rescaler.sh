#!/bin/bash

# Source common wrapper functions
COMMONDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../common" && pwd )"
source $COMMONDIR/common_wrapper.sh
source $COMMONDIR/common_inputs.sh

# Command line
while [ $# -gt 0 ]; do
    case $1 in
        -site|--site) site=$2; verify_not_null $site; shift 2;;
        -param*|--param*) params=$2; verify_not_null $params; shift 2;;
        *) echo "Do not recognize command line option: $1" 1>&2; usage;;
    esac
done

if [ -z "$site" ] || [ ! -f "$params" ]; then
    usage
fi

params=$( readlink -f $params )
utils=$( readlink -f ../../utils )
PATH=$PATH:$utils:$PWD
aggr_mask_directory=$( get_param aggr_mask_directory )
level=$( get_param rescaler:level )
mkfile=${aggr_mask_directory}/$level.mask.nc4

swift -tc.file tc.data -sites.file ${site}.xml -config swift.properties rescaler.swift \
      -mkfile=$mkfile \
      -agglvl=$level \
      -params=$params
cleanup $?
