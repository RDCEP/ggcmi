#!/bin/bash

export PYTHONPATH=$PYTHONPATH:/project/joshuaelliott/ggcmi/utils
rescaler.py "$@" 2>&1
exit 0
