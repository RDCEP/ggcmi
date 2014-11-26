#!/bin/bash

root=/project/ggcmi

# Header
echo indir metricsdir agglvl outdir

# gadm0
echo $root/AgMIP.output/processed/biascorr/gadm0/faostat $root/AgMIP.output/processed/multimetrics/gadm0/faostat gadm0 $root/AgMIP.output/processed/modelensemble/gadm0/faostat

# fpu
echo $root/AgMIP.output/processed/biascorr/fpu/ray $root/AgMIP.output/processed/multimetrics/fpu/ray fpu $root/AgMIP.output/processed/modelensemble/fpu/ray
echo $root/AgMIP.output/processed/biascorr/fpu/iizumi $root/AgMIP.output/processed/multimetrics/fpu/iizumi fpu $root/AgMIP.output/processed/modelensemble/fpu/iizumi

# kg
echo $root/AgMIP.output/processed/biascorr/kg/ray $root/AgMIP.output/processed/multimetrics/kg/ray kg $root/AgMIP.output/processed/modelensemble/kg/ray
echo $root/AgMIP.output/processed/biascorr/kg/iizumi $root/AgMIP.output/processed/multimetrics/kg/iizumi kg $root/AgMIP.output/processed/modelensemble/kg/iizumi

# global
echo $root/AgMIP.output/processed/biascorr/global/faostat $root/AgMIP.output/processed/multimetrics/global/faostat global $root/AgMIP.output/processed/modelensemble/global/faostat
