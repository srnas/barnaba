from __future__ import absolute_import, division, print_function

import barnaba.commandline as cl

import barnaba.enm as enm
import os
from comp_mine import comp
import numpy as np

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd

fname_prot= "%s/test/data/2n82.pdb" % cwd

def test_enm_0():
    """Tests command line interface"""
    cl.main(['ENM','--pdb',fname,'-o',"%s/enm_00.test" % outdir,'--ntop','4','--shape'])

    comp("%s/enm_00.test.ENM.eigvecs.out" % refdir)
    comp("%s/enm_00.test.ENM.eigvals.out" % refdir)
    comp("%s/enm_00.test.ENM.shape.out" % refdir)

def test_enm_1():
    '''Tests C2-C2 calculation'''
    
    # initialize class. Only PDB are accepted
    SBP_enm = enm.Enm(fname,sele_atoms=["C2","C1\'","P"])
    
    # print eigenvectors
    evecs = SBP_enm.print_evec(3)
    fh = open("%s/enm_01.test.dat" % outdir,'w')    
    fh.write(evecs)
    fh.close()

    # print eigenvalues 
    evals = SBP_enm.print_eval()
    fh = open("%s/enm_02.test.dat" % outdir,'w')    
    fh.write(evals)
    fh.close()

    # print mean square fluctuations
    fluc = SBP_enm.get_MSF()
    np.savetxt("%s/enm_03b.test.dat" % outdir,fluc,fmt='%.6e')
    
    # print C2'-C2' fluctuations
    fluc,res =  SBP_enm.c2_fluctuations()
    fh = open("%s/enm_03.test.dat" % outdir,'w')    
    stri = "# %19s %s \n" % ("Residues","Fluctuations")
    for i in range(len(fluc)):
        stri +=  "%10s/%-10s %.6e \n" % (res[i],res[i+1],fluc[i])
    fh.write(stri)
    fh.close()

    comp("%s/enm_01.test.dat" % refdir)
    comp("%s/enm_02.test.dat" % refdir)
    comp("%s/enm_03.test.dat" % refdir)
    comp("%s/enm_03b.test.dat" % refdir)

def test_enm_2():
    '''Tests different bead choices'''
    
    # initialize class. Only PDB are accepted
    S_enm = enm.Enm(fname,sele_atoms=["C1\'"],cutoff=1.5)
    SB_enm = enm.Enm(fname,sele_atoms=["C2","C1\'"],cutoff=1.1)
    
    # print eigenvalues S-ENM 
    evals = S_enm.print_eval()
    fh = open("%s/enm_04.test.dat" % outdir,'w')    
    fh.write(evals)
    fh.close()

    # print eigenvalues SB-ENM
    evals = SB_enm.print_eval()
    fh = open("%s/enm_05.test.dat" % outdir,'w')    
    fh.write(evals)
    fh.close()
    
    comp("%s/enm_04.test.dat" % refdir)
    comp("%s/enm_05.test.dat" % refdir)

def test_enm_3():
    '''Tests sparse diagonalization'''
    
    # initialize class. Only PDB are accepted
    AA_enm = enm.Enm(fname,sele_atoms="AA",sparse=True,cutoff=0.7)
    
    # print eigenvalues 
    evals = AA_enm.print_eval()
    fh = open("%s/enm_06.test.dat" % outdir,'w')    
    fh.write(evals)
    fh.close()

    # print eigenvectors
    evecs = AA_enm.print_evec(3)
    fh = open("%s/enm_07.test.dat" % outdir,'w')    
    fh.write(evecs)
    fh.close()

    # print C2'-C2' fluctuations
    fluc,res =  AA_enm.c2_fluctuations()
    fh = open("%s/enm_08.test.dat" % outdir,'w')    
    stri = "# %19s %s \n" % ("Residues","Fluctuations")
    for i in range(len(fluc)):
        stri +=  "%10s/%-10s %.6e \n" % (res[i],res[i+1],fluc[i])
    fh.write(stri)
    fh.close()
    
    comp("%s/enm_06.test.dat" % refdir)
    comp("%s/enm_07.test.dat" % refdir)
    comp("%s/enm_08.test.dat" % refdir)

    

def test_enm_4():
    '''Tests RNA+protein complex'''

    # initialize class. Only PDB are accepted                                                                                                                                             
    SBP_enm = enm.Enm(fname_prot,sparse=False)

    # print eigenvalues
    evals = SBP_enm.print_eval()
    fh = open("%s/enm_09.test.dat" % outdir,'w')
    fh.write(evals)
    fh.close()

    # print eigenvectors
    evecs = SBP_enm.print_evec(3)
    fh = open("%s/enm_10.test.dat" % outdir,'w')
    fh.write(evecs)
    fh.close()

    # print C2'-C2' fluctuations
    fluc,res =  SBP_enm.c2_fluctuations()
    fh = open("%s/enm_11.test.dat" % outdir,'w')
    stri = "# %19s %s \n" % ("Residues","Fluctuations")
    for i in range(len(fluc)):
        stri +=  "%10s/%-10s %.6e \n" % (res[i],res[i+1],fluc[i])
    fh.write(stri)
    fh.close()
    comp("%s/enm_09.test.dat" % refdir)
    comp("%s/enm_10.test.dat" % refdir)
    comp("%s/enm_11.test.dat" % refdir)

def test_enm_5():
    '''Tests RNA+protein complex (sparse)'''

    # initialize class. Only PDB are accepted
    AA_enm = enm.Enm(fname_prot,sele_atoms="AA",sparse=True,cutoff=0.7)

    # print eigenvalues
    evals = AA_enm.print_eval()
    fh = open("%s/enm_12.test.dat" % outdir,'w')
    fh.write(evals)
    fh.close()

    # print eigenvectors
    evecs = AA_enm.print_evec(3)
    fh = open("%s/enm_13.test.dat" % outdir,'w')
    fh.write(evecs)
    fh.close()


    comp("%s/enm_12.test.dat" % refdir)
    comp("%s/enm_13.test.dat" % refdir)

    
def test_enm_6():
    """Tests modes output"""
    # initialize class. Only PDB are accepted                                                                                                                                                          
    #SBP_enm = enm.Enm(fname,sparse=True)
    SBP_enm = enm.Enm(fname,sparse=True)

    import mdtraj as md
    traj_mode=SBP_enm.get_mode_traj(6)
    traj_mode.save_pdb("%s/enm_14.test.pdb" % outdir)
    
    comp("%s/enm_14.test.pdb" % refdir)

