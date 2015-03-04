#!/bin/bash

set -e

# Search for GNRA tetraloop in 1S72 PDB file. -f option accepts multiple files and/or concatenated PDB
# 1. The maximum ERMSD treshold value can be specified using the --treshold option. As a rule of thumb, 
#    values between 0.5-1.0 should work in most cases.  
# 2. The maximum number of allowed bulged bases (max 2) can be specified using the --bulge option 

../baRNAba --name 3_EXAMPLE_1 SS_MOTIF --query DATA/GNRA.pdb -f DATA/1S72.pdb --bulge 1 --treshold 0.6
