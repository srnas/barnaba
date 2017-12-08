
Introduction
============

BaRNAba is a tool for analyzing RNA three-dimensional structures.
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
To perform stop motion modeling analysis (SMM), pyemma is required too.

We highly recommend to install the above packages using conda (http://conda.pydata.org/docs/index.html)

  conda install -c omnia mdtraj
  conda install -c omnia pyemma

Installation
-------------
You can obtain barnaba using git:

    git clone -b library git://github.com/srnas/barnaba.git

or download a zip file from the web:

   https://github.com/srnas/barnaba/tree/library

Usage
------------
BaRNAba can be either used as a python package or using the command-line interface.
Here's some minimal examples. You can find more elaborate examples at XXX
0.  minimal help
    ./baRNAba --help 

1. Calculate the ERMSD between structures
   ./baRNAba ERMSD --ref sample1.pdb --pdb sample2.pdb sample3.pdb ...
   
   trajectories can be provided as well, by specifying a topology file
   ./baRNAba ERMSD --ref sample1.pdb --top sample1.pdb --trj samples.xtc 

2. Calculate ESCORE
   ./baRNAba ESCORE --ff 1S72.pdb --pdb sample1.pdb ...

3. Find hairpin loop motif
   ./baRNAba SS_MOTIF --query motif.pdb --pdb file1.pdb file2.pdb 

4. Find double stranded motif. l1 and l2 are the lengths of the two strands
   ./baRNAba DS_MOTIF --query motif.pdb --pdb file1.pdb file2.pdb ... --l1 l1 --l2 l2

5. Annotate structures/trajectories according to the Leontis/Westhof classification.  
   ./baRNAba ANNOTATE --pdb file1.pdb file2.pdb ...

6. Calculate dihedral backbone angles and sugar pucker
   ./baRNAba TORSION --pdb file1.pdb file2.pdb --hread 

7. Calculate J-couplings (H1'H2', H2'H3' H3'H4', H4'H5', H4'H5'',1H5P,2H5P,H3P_+)
   ./baRNAba JOCUPLING --pdb file1.pdb file2.pdb ... --jcoupling --hread 

8. Calculate NOE signal
   ./baRNAba NOE --pdb file1.pdb 

9. Calculate elastic network models for RNA and predict SHAPE reactivity 
   ./baRNAba ENM --pdb file1.pdb 












