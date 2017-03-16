#!/bin/bash

export PATH=$PATH:$PWD

site=$1
if [ -z "$site" ]; then
    site="sandyb"
fi

swift -sites.file ${site}.xml -tc.file tc.data aggregator.swift

# Remove run directories if Swift finishes with no errors
if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
