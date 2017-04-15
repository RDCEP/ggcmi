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

swift -tc.file tc.data -sites.file ${site}.xml -config swift.properties rescaler.ensemble.swift -params=$params
cleanup $?
