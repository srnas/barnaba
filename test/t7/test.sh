#!/bin/bash

set -e

# calculate TORSION and pucker angles.

$BARNABA --name 7_EXAMPLE_TOR TORSION  -f $DATA/1S72.pdb
$BARNABA --name 7_EXAMPLE_PUC TORSION  -f $DATA/1S72.pdb



