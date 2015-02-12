#!/bin/bash

# only process gadm0, faostat, fixed_mirca_mask for now

edir=/project/ggcmi/AgMIP.output/processed/modelensemble/gadm0/faostat/fixed_mirca_mask

indir=/project/ggcmi/AgMIP.output/processed/rescaled/gadm0/faostat/fixed_mirca_mask

mkfile=/project/ggcmi/AgMIP.output/processed/masks/aggr/gadm0.mask.nc4

crmthd=2,0,2 # dt, mp, cr

agglvl=gadm0

outdir=/project/ggcmi/AgMIP.output/processed/rescaled.ensemble/gadm0/faostat/fixed_mirca_mask

# Header
echo infile indir mkfile crmthd agglvl outfile

mkdir -p $outdir

for f in `ls $edir`; do
   echo $edir/$f $indir $mkfile $crmthd $agglvl $outdir/${f/ensemble/rescaled}
done
