#!/bin/bash

set -e

# Search for SARCIN motif in 1S72 PDB file. -f option accepts multiple files and/or concatenated PDB
# the maximum number of allowed per-strand inserted can be specified using the -bulge option 
$BARNABA DS_MOTIF --query $DATA/SARCIN.pdb --pdb $DATA/1S72.pdb --l1 8 --l2 7 --bulge 1 --treshold 0.65 --dump -o DSMOTIF
