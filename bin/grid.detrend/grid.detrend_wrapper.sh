#!/bin/bash

INFILE=$1
LATIDX=$2
LONIDX=$3
NUMLAT=$4
NUMLON=$5
NBTLAT=$6
NBTLON=$7
VRLIST=$8
OUTDIR=$9

mkdir -p $OUTDIR/$LATIDX

SLAT=$(( NUMLAT / NBTLAT ))
SLON=$(( NUMLON / NBTLON ))

LATIDX1=$(( (LATIDX - 1) * SLAT ))
if [ $LATIDX = $NBTLAT ]; then
   LATIDX2=$(( NUMLAT - 1 ))
else
   LATIDX2=$(( LATIDX * SLAT - 1 ))
fi

LONIDX1=$(( (LONIDX - 1) * SLON ))
if [ $LONIDX = $NBTLON ]; then
   LONIDX2=$(( NUMLON - 1 ))
else
   LONIDX2=$(( LONIDX * SLON - 1 ))
fi

TMPFILE=tmp_${LATIDX}_${LONIDX}.nc4

ncks -h -d lat,$LATIDX1,$LATIDX2 -d lon,$LONIDX1,$LONIDX2 $INFILE $TMPFILE

OUTFILE=$OUTDIR/$LATIDX/$(basename $INFILE.$LONIDX)

/project/joshuaelliott/ggcmi/bin/grid.detrend/grid.detrend.py -i $TMPFILE -v $VRLIST -o $OUTFILE

rm $TMPFILE
