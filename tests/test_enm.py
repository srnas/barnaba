
import barnaba.enm as enm

fname = "data/sample1.pdb"

# initialize class. Only PDB are accepted
SBP_enm = enm.Enm(fname,sele_atoms=["C2","C1\'","P"])

# print eigenvectors
evecs = SBP_enm.print_evec(3)
fh = open("enm_evec.test.dat",'w')
fh.write(evecs)
fh.close()

# print eigenvalues 
evals = SBP_enm.print_eval()
fh = open("enm_eval.test.dat",'w')
fh.write(evals)
fh.close()

# print C2'-C2' fluctuations
fluc,res =  SBP_enm.c2_fluctuations()
fh = open("enm_fluctuations.test.dat",'w')
stri = "# %19s %s \n" % ("Residues","Fluctuations")
for i in range(len(fluc)):
    stri +=  "%10s/%-10s %.6e \n" % (res[i],res[i+1],fluc[i])
fh.write(stri)
fh.close()
