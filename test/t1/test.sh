#!/bin/bash

set -e

# Calculate ERMSD distance between sample1.pdb and sample2.pdb. -f option accepts multiple files
$BARNABA --name ERMSD_EXAMPLE_1 ERMSD --pdb $DATA/sample1.pdb -f $DATA/sample2.pdb

# Calculate ERMSD distance between sample1.pdb and models in samples.pdb
$BARNABA --name ERMSD_EXAMPLE_2 ERMSD --pdb $DATA/sample1.pdb -f $DATA/samples.pdb

# Calculate ERMSD distance between sample1.pdb and sample2.pdb. 
# Additionally calculate per-residue ERMSD flucutations
$BARNABA --name ERMSD_EXAMPLE_3 ERMSD --pdb $DATA/sample1.pdb -f $DATA/sample2.pdb --ermsf


