#!/bin/bash

utils=$( readlink -f ../../utils )
PATH=$PATH:$utils:$PWD

site=$1
if [ -z "$site" ]; then
    site="sandyb"
fi

swift -tc.file tc.data -sites.file ${site}.xml -config swift.properties multimetrics.ensemble.swift
rc=$?
echo Cleaning up, please wait
sleep 5

if [ "$rc" -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
