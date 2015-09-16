#!/bin/bash

set -e

# Calculate ESCORE for all models in samples.pdb using pdb 1S72 as "force field".
# -f option accepts multiple files
../baRNAba --name ESCORE_EXAMPLE_1 ESCORE --force-field DATA/1S72.pdb -f DATA/samples.pdb


