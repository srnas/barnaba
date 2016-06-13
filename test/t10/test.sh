#!/bin/bash

set -e

# calculate backbone and pucker angles
$BARNABA SMM --gvec w_repetition.gvec --ref ref.gvec --ncluster 100 -o plain
awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' plain.SMM.rna.smm.dat > plain_smm.dat

$BARNABA SMM --gvec no_repetition.gvec --weight weights.dat --ref ref.gvec --ncluster 100 -o weight
awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' plain.SMM.rna.smm.dat > weight_smm.dat



