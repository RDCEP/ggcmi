#!/bin/bash

PYTHONPATH=$PYTHONPATH:/project/joshuaelliott/ggcmi2/utils
aggregator.py "$@" 2>&1
exit 0
