#!/bin/bash

set -e

# Calculate ENM as described in 
# Pinamonti et al "Elastic Network Model for RNA: a comparative assessment with MD and SHAPE experiments" (2015) 

$BARNABA  ENM --pdb $DATA/sample1.pdb --type SBP -o ENM_1
$BARNABA  ENM --pdb $DATA/sample1.pdb --type AA --cutoff 0.7 -o ENM_2

#$BARNABA  ENM --pdb $DATA/sample1.pdb --type SBP -o ENM_2


