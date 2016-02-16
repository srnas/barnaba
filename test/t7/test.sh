#!/bin/bash

set -e

# calculate backbone and pucker angles
$BARNABA  TORSION  --pdb $DATA/1S72.pdb --hread -o T1
$BARNABA  TORSION  --top $DATA/sample1.pdb --trj $DATA/samples.xtc --hread -o T2

# calculate j couplings
$BARNABA  JCOUPLING --top $DATA/sample1.pdb --trj $DATA/samples.xtc --hread -o T3





