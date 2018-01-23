.. image:: https://travis-ci.org/srnas/barnaba.svg?branch=regression-tests
    :target: https://travis-ci.org/srnas/barnaba

Introduction
============

BaRNAba is a tool for analyzing RNA three-dimensional structures and simulations. It supports several formats including pdb, xtc, trr, dcd, binpos, netcdf, mdcrd, prmtop, and more, as it uses the MDtraj library.
BaRNAba has been developed by Sandro Bottaro with the crucial help of Giovanni Bussi and Giovanni Pinamonti.
This is what you can do with baRNAba:

1. Calculate eRMSD for nucleic acids
2. Calculate RMSD and do RMSD superpositions
3. Search for single stranded RNA motifs/fragments in the PDB database or in simulations
4. Search for double stranded RNA motifs/fragments in the PDB database or in simulations
5. Annotate PDB structures and trajectories with the Leontis-Westhof classification
6. Cluster nucleic acids structures using the eRMSD as a metric distance
7. Calculate elastic network models for nucleic acids and nucleic acids/protein complexes (Pinamonti et al. Nucleic Acids Research, 2015)
8. Perform a Stop-Motion analysis (Bottaro et al. Nucleic Acids Research, 2016)
9. Calculate backbone and pucker torsion angles in a PDB structure or trajectory
10. Back-calculate 3J scalar couplings from PDB structure or trajectory
11. Score three-dimensional structures using eSCORE (Bottaro et al, Nucleic Acids Research)

For bugs, questions or comments contact Sandro at sandro dot bottaro (guess) gmail dot com

If you use baRNAba in your work,  please cite the following paper:

S.Bottaro, F. di Palma, G.Bussi. The role of nucleobases 
in RNA structure and dynamics.  Nucleic Acids Research (2014)

Requirements
-------------
baRNAba requires:
   - Python 2.7.x
   - Numpy
   - Scipy
   - Mdtraj
     
baRNAba requires mdtraj (http://mdtraj.org/) for manipulating structures and trajectories
To perform cluster analysis, scikit-learn is required too.

MDtraj can be installed using pip:

    pip install mdtraj

Installation
-------------
You can obtain barnaba using git:

    git clone git://github.com/srnas/barnaba.git

or download a zip file from the web:

   https://github.com/srnas/barnaba/tree/library

then move to barnaba directory and run the command

   pip install -e .


Usage
------------
BaRNAba can be either used as a Python library or as a commandline tool.
A number of Notebook examples can be found in the examples directory. Alternatively, the command-line interface can be found in the bin directory. Here's a minimal how-to

0.  minimal help  

    ./baRNAba.py --help

1. Calculate the ERMSD between structures  

   ./baRNAba.py ERMSD --ref ../test/data/sample1.pdb --pdb ../test/data/sample2.pdb
  
trajectories can be provided as well, by specifying a topology file  

   ./baRNAba.py ERMSD --ref ../test/data/sample1.pdb --top ../test/data/sample1.pdb --trj ../test/data/samples.xtc  

other accepted options are shown in a function-specific help  

   ./baRNAba.py ERMSD --help
  
2. Calculate the RMSD between structures  
  
   ./baRNAba.py RMSD --ref ../test/data/sample1.pdb --pdb ../test/data/sample2.pdb --dump
   
3. Find single stranded motif  
  
   ./baRNAba.py SS_MOTIF --query ../test/data/GNRA.pdb --pdb ../test/data/1S72.pdb   
   
4. Find double stranded motif. l1 and l2 are the lengths of the two strands
  
   ./baRNAba.py DS_MOTIF --query ../test/data/SARCIN.pdb --pdb ../test/data/1S72.pdb --l1 8 --l2 7  
 
5. Annotate structures/trajectories according to the Leontis/Westhof classification.
   
   ./baRNAba.py ANNOTATE --pdb ../test/data/SARCIN.pdb  

6. Calculate backbone/sugar/pseudorotation angles
    
   ./baRNAba.py TORSION --pdb ../test/data/GNRA.pdb --backbone --sugar --pucker 
 

7. Calculate J-couplings 

   ./baRNAba.py JCOUPLING --pdb ../test/data/sample1.pdb 

9. Calculate elastic network models for RNA and predict SHAPE reactivity. NB: only works with PDB.
   
   ./baRNAba.py ENM --pdb ../test/data/GNRA.pdb --shape

10. Calculate relative positions between bases R_ij  ang G vectors for pairs within ellipsoidal cutoff  

   ./baRNAba.py DUMP --pdb ../test/data/GNRA.pdb --dumpG --dumpR  

11. Extract fragments from structures with a given sequence. NB: only works with PDB.  

    ./baRNAba.py SNIPPET --pdb ../test/data/1S72.pdb  --seq NNGNRANN
 
12. Calculate ESCORE  
    ./baRNAba.py ESCORE --ff ../test/data/1S72.pdb --pdb ../test/data/sample1.pdb













