#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file apps -sites.file midway.xml -config swift.properties biascorrect.swift

if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
