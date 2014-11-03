#!/bin/bash

root=/project/joshuaelliott/ggcmi

# Header
echo indir reffile agglvl outdir

# gadm0
echo $root/processed/biascorr/gadm0/faostat $root/reference/faostat/faostat.1961-2012.gadm0.nc4 gadm0 $root/processed/multimetrics/gadm0/faostat

# fpu
echo $root/processed/biascorr/fpu/ray $root/reference/ray/ray.1961-2008.fpu.nc4 fpu $root/processed/multimetrics/fpu/ray
echo $root/processed/biascorr/fpu/iizumi $root/reference/iizumi/iizumi.1982-2006.fpu.nc4 fpu $root/processed/multimetrics/fpu/iizumi

# kg
echo $root/processed/biascorr/kg/ray $root/reference/ray/ray.1961-2008.kg.nc4 kg $root/processed/multimetrics/kg/ray
echo $root/processed/biascorr/kg/iizumi $root/reference/iizumi/iizumi.1982-2006.kg.nc4 kg $root/processed/multimetrics/kg/iizumi
