#!/bin/bash

root=/project/joshuaelliott/ggcmi

# Header
echo indir reffile agglvl outdir

# gadm0
echo $root/processed/aggs/gadm0 $root/reference/faostat/faostat.1961-2012.gadm0.nc4 gadm0 $root/processed/biascorr/gadm0/faostat

# fpu
echo $root/processed/aggs/fpu $root/reference/ray/ray.1961-2008.fpu.nc4 fpu $root/processed/biascorr/fpu/ray
echo $root/processed/aggs/fpu $root/reference/iizumi/iizumi.1982-2006.fpu.nc4 fpu $root/processed/biascorr/fpu/iizumi

# kg
echo $root/processed/aggs/kg $root/reference/ray/ray.1961-2008.kg.nc4 kg $root/processed/biascorr/kg/ray
echo $root/processed/aggs/kg $root/reference/iizumi/iizumi.1982-2006.kg.nc4 kg $root/processed/biascorr/kg/iizumi
