
import barnaba.enm as enm
import os
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd

def test_enm_1():
    
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
