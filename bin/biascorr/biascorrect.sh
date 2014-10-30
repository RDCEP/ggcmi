#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file apps -sites.file midway.xml -config swift.properties biascorrect.swift

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
