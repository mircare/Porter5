#!/bin/bash

# Transforms a HHblits psi file into a flat MA file containing in the first line
# the number of alignments and the sequences in the following lines


ff=`echo $1 | sed -e 's/\.psi/\.flatpsi/g'`
wc -l $1 |awk '{print $1}' > $ff
awk '{print $2}' $1 >> $ff
sed -i'' -e 's/-/./g' $ff
sed -i'' -e '2s/\.//g' $ff
