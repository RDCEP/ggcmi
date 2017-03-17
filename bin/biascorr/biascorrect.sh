#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

site=$1
if [ -z "$site" ]; then
    site="sandyb"
fi

swift -sites.file ${site}.xml -tc.file tc.data -config swift.properties biascorrect.swift

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
