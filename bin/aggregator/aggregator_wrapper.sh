#!/bin/bash

PATH=$PATH:/project/joshuaelliott/ggcmi/utils
aggregator.py "$@" 2>&1 
exit 0
