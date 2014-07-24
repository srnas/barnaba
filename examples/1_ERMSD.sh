#!/bin/bash

set -e

# Calculate ERMSD distance between sample1.pdb and sample2.pdb. -f option accepts multiple files
baRNAba --name ERMSD_EXAMPLE_1 ERMSD --pdb DATA/sample1.pdb -f DATA/sample2.pdb

# Calculate ERMSD distance between sample1.pdb and models in samples.pdb
baRNAba --name ERMSD_EXAMPLE_2 ERMSD --pdb DATA/sample1.pdb -f DATA/samples.pdb


# Calculate ERMSD distance between sample1.pdb and samples in gromacs trajectory.
# XTC and GRO formats are supported. 
# when analyzing trajectories, trajconv executable and reference PDB file must be specified!
baRNAba --trjconv trjconv --name ERMSD_EXAMPLE_3 ERMSD --pdb DATA/sample1.pdb -f DATA/trj.xtc


# Calculate ERMSD distance between sample1.pdb and sample2.pdb. 
# Additionally calculate per-residue ERMSD flucutations
baRNAba --name ERMSD_EXAMPLE_4 ERMSD --pdb DATA/sample1.pdb -f DATA/sample2.pdb --ermsf

# Calculate ERMSD distance between sample1.pdb and sample2.pdb. 
# Use contact-map-like ERMSD 
baRNAba --name ERMSD_EXAMPLE_5 ERMSD --pdb DATA/sample1.pdb -f DATA/sample2.pdb --type scalar --cutoff 2.0

