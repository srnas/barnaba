#!/bin/bash

set -e

# calculate backbone and pucker angles

$BARNABA NOE --top $DATA/sample1.pdb --trj $DATA/samples.xtc -o NOE --nbins 5




