#!/bin/bash

set -e

# Calculate ESCORE for all models in samples.pdb using pdb 1S72 as "statistical force field".

$BARNABA ESCORE --ff $DATA/1S72.pdb --top $DATA/sample1.pdb --trj $DATA/samples.xtc -o ESCORE_1

$BARNABA ESCORE --ff $DATA/1S72.pdb --pdb $DATA/sample1.pdb -o ESCORE_2


