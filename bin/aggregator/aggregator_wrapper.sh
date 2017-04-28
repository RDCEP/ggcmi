#!/bin/bash

PYTHONPATH=$PYTHONPATH:/project/joshuaelliott/ggcmi/utils
aggregator.py "$@" 2>&1
exit 0
