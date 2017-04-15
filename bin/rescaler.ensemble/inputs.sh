#!/bin/bash

# Parameter file
params=$1
if [ -z "$params" ]; then
    echo "Usage: $0 <params>"
    exit 1
fi

COMMONDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../common" && pwd )"
source $COMMONDIR/common_inputs.sh

weight_directory=$( get_param weight_directory )
rescaled_directory=$( get_param rescaled_directory )
area_long=$( area_to_long $( get_param rescaler:area ))
modelensemble_directory=$( get_param modelensemble_directory )
reference=$( get_param rescaler:reference )
level=$( get_param rescaler:level )
rescaledensemble_directory=$( get_param rescaledensemble_directory )
aggr_mask_directory=$( get_param aggr_mask_directory)

# only process gadm0, faostat, fixed_mirca_mask for now
edir=$modelensemble_directory/$level/$reference/$area_long
indir=$rescaled_directory/$level/$reference/$area_long
mkfile=$aggr_mask_directory/$level.mask.nc4
level=$( get_param rescaler:level )
outdir=$rescaledensemble_directory/$level/$reference/$area_long

# Header
echo infile indir mkfile agglvl outfile

mkdir -p $outdir

for f in `ls $edir`; do
   echo $edir/$f $indir $mkfile $level $outdir/${f/ensemble/rescaled}
done
