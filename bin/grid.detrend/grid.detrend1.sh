#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils

swift -tc.file tc.data -sites.file midway.xml -config swift.properties grid.detrend1.swift

if [ $? -eq 0 ]; then
   rm -rf run???
   rm finder.out
fi
