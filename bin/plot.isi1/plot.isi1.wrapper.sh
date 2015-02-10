#!/bin/bash

export PATH=$PATH:$PWD

swift -sites.file midway.xml -tc.file tc.data plot.isi1.swift \
      -plots=blmap,dymap,box                                  \
      -crops=maize,wheat,soy,rice,all

# Remove run directories if Swift finishes with no errors
if [ $? -eq 0 ]; then
   echo Removing run directory . . .
   rm -rf run???
fi
