#!/bin/bash

INDIR=$1
LATIDX=$2
INFILE=$3

for f in $INDIR/$LATIDX/$INFILE.*; do
   ncpdq -O -h -a lon,time $f $f
   ncks -O -h --mk_rec_dim lon $f $f
done

OUTFILE=$INDIR/$INFILE.$LATIDX

ncrcat -h $INDIR/$LATIDX/$INFILE.* $OUTFILE

ncpdq -O -h -a lat,lon $OUTFILE $OUTFILE
ncks -O -h --mk_rec_dim lat $OUTFILE $OUTFILE
