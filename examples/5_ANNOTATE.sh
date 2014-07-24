#!/bin/bash

set -e

# Annotate structures using the Leontis-Westhof classification. 
# results are similar (but not identical) to MC annotate. 

baRNAba --name 5_EXAMPLE ANNOTATE -f DATA/1S72.pdb
