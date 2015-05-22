#!/bin/bash

set -e

# Calculate ENM as described in 
# Pinamonti et al "Elastic Network Model for RNA: a comparative assessment with MD and SHAPE experiments" (2015) 

$BARNABA --name 6_EXAMPLE_SBP ENM -f $DATA/sample1.pdb --type SBP

$BARNABA --name 6_EXAMPLE_AA ENM -f $DATA/sample1.pdb --type AA --cutoff 7.0
