#!/bin/bash

root=/project/joshuaelliott/ggcmi

# Header
echo indir metricsdir agglvl outdir

# gadm0
echo $root/processed/biascorr/gadm0/faostat $root/processed/multimetrics/gadm0/faostat gadm0 $root/processed/modelensemble/gadm0/faostat

# fpu
echo $root/processed/biascorr/fpu/ray $root/processed/multimetrics/fpu/ray fpu $root/processed/modelensemble/fpu/ray
echo $root/processed/biascorr/fpu/iizumi $root/processed/multimetrics/fpu/iizumi fpu $root/processed/modelensemble/fpu/iizumi

# kg
echo $root/processed/biascorr/kg/ray $root/processed/multimetrics/kg/ray kg $root/processed/modelensemble/kg/ray
echo $root/processed/biascorr/kg/iizumi $root/processed/multimetrics/kg/iizumi kg $root/processed/modelensemble/kg/iizumi
