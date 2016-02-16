#!/bin/bash

set -e

# Annotate structures using the Leontis-Westhof classification. 
# The annotation is not accurate because  - especially for stacking - it is not well defined.
# if you want to get accurate annotation with hydrogen bonds, etc. I highly reccommend to use X3DNA.

$BARNABA ANNOTATE --pdb $DATA/1S72.pdb --hread -o ANNO_1
$BARNABA ANNOTATE --pdb $DATA/1tra.pdb --hread -o ANNO_2 --pymol
