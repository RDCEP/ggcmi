#!/bin/bash

# =================================
# GGCMI PHASE 1 PROCESSING PIPELINE
# =================================

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

# Downscaled
cd /project/joshuaelliott/ggcmi/bin/downscaler
./downscaler.sh

# Rescaled ensemble
cd /project/joshuaelliott/ggcmi/bin/rescaler.ensemble
./rescaler.ensemble.sh

cd $dir
