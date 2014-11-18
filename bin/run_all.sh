#!/bin/bash

dir=$PWD

# Aggregator
cd /project/joshuaelliott/ggcmi/bin/aggregator
/project/joshuaelliott/ggcmi/bin/aggregator/aggregator.sh

# Bias correcter
cd /project/joshuaelliott/ggcmi/bin/biascorr
/project/joshuaelliott/ggcmi/bin/biascorr/biascorrect.sh

cd $dir
