#!/bin/bash

for c in maize rice soy wheat; do
   for d in 2010 2020 2030 2040 2050; do
      for co in co2 noco2; do
         echo $c, $d, $co . . .
         ./relyieldcsv.py -i /project/ggcmi/isi1/processed/isi1.long.agg -c $c --co2 $co -d $d -o .
      done
   done
done
