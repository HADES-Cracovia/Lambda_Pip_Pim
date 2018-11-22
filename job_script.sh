#!/bin/bash

date

#. /lustre/nyx/hades/user/rlalik/hades/pp45/profile.sh
. /lustre/nyx/hades/user/knowakow/HYPERON/Lambda1520/Lambda1520_ic/profile.sh

echo file=$pattern
echo events=$events
echo odir=$odir

root -b -q

cd /lustre/nyx/hades/user/knowakow/HYPERON/Lambda1520/Lambda1520_ic

date

ofile=$(basename $pattern .root)_ana.root
fname=$(echo $ofile | cut -d '_' -f 3)

[ ! -d  $odir/$fname ] && mkdir $odir/$fname

echo ./anaLambda1520 -e -1 $pattern -o $odir/$fname/$ofile
time ./anaLambda1520 -e -1 $pattern -o $odir/$fname/$ofile
