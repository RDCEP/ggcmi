#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

site=$1
if [ -z "$site" ]; then
    site="sandyb"
fi

swift -tc.file tc.data -sites.file ${site}.xml -config swift.properties multimetrics.swift

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
