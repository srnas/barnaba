#!/bin/bash

set -e

# Calculate ERMSD distance between sample1.pdb and sample2.pdb. -f option accepts multiple files
$BARNABA ERMSD --ref $DATA/sample1.pdb --pdb $DATA/sample2.pdb --per-res -o ERMSD_1

# Calculate ERMSD distance between sample1.pdb and samples.xtc trajectory
$BARNABA ERMSD --ref $DATA/sample1.pdb --trj $DATA/samples.xtc --top  $DATA/sample1.pdb  -o ERMSD_2



