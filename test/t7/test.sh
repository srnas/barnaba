#!/bin/bash

set -e

# calculate backbone and pucker angles

$BARNABA --name 7_EXAMPLE_TOR TORSION  -f $DATA/1S72.pdb --pucker --hread
$BARNABA --name 7_EXAMPLE_TOR TORSION  -f $DATA/1S72.pdb --backbone --hread

# calculate j couplings 
$BARNABA --name 7_EXAMPLE_TOR TORSION -f $DATA/samples.pdb --jcoupling --hread




