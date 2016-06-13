#!/bin/bash

set -e

# calculate backbone and pucker angles
$BARNABA SMM --gvec w_repetition.gvec --ref ref.gvec --ncluster 100 -o plain
$BARNABA SMM --gvec no_repetition.gvec --weight weights.dat --ref ref.gvec --ncluster 100 -o weight




