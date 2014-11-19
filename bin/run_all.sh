#!/bin/bash

dir=$PWD

# Aggregator
cd /project/joshuaelliott/ggcmi/bin/aggregator
./aggregator.sh

# Bias correcter
cd /project/joshuaelliott/ggcmi/bin/biascorr
./biascorrect.sh

# Multimetrics
cd /project/joshuaelliott/ggcmi/bin/multimetrics
./multimetrics.sh

# Model ensemble
cd /project/joshuaelliott/ggcmi/bin/modelensemble
./modelensemble.sh

# Ensemble multimetrics
cd /project/joshuaelliott/ggcmi/bin/multimetrics.ensemble
./multimetrics.ensemble.sh

# Rescaled
cd /project/joshuaelliott/ggcmi/bin/rescaler
./rescaler.sh

cd $dir
