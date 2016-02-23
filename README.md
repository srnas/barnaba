[![Build Status](https://travis-ci.org/srnas/barnaba.svg)](https://travis-ci.org/srnas/barnaba)

###### INTRODUCTION ##############3
BaRNAba is a tool for analyzing RNA three-dimensional structures.
BaRNAba has been developed by Sandro Bottaro with the crucial help of Giovanni Bussi.
For bugs, questions or comments contact Sandro at sbottaro@sissa.it. 
If you use baRNAba in your work,  please cite the following paper:

S.Bottaro, F. di Palma, G.Bussi. The role of nucleobases 
in RNA structure and dynamics.  Nucleic Acids Research (2014)

###### REQUIREMENTS ##############
baRNAba requires mdtraj (http://mdtraj.org/) for manipulating structures and trajectories
To perform stop motion modeling analysis (SMM), pyemma is required too. 
We highly recommend to install the above packages using conda (http://conda.pydata.org/docs/index.html)

  conda install -c omnia mdtraj
  conda install -c omnia pyemma
  
###### INSTALLATION ###############
You can obtain barnaba using git:

    git clone git://github.com/srnas/barnaba.git

or download a zip file from the web:

   https://github.com/srnas/barnaba


### USAGE

* minimal help
  ./baRNAba --help 

Currently, baRNAba can perform different tasks:

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












