
import barnaba.enm as enm

fname = "data/sample1.pdb"

SBP_enm = enm.Enm(fname,sele_atoms=["C2","C1\'","P"])

print SBP_enm.print_evec(3)

exit()
fh = open("couplings.dat",'w')
stri = ""
for e in range(len(rr)):
    stri += "%20s " % (rr[e])
    for k in range(10):
        stri += " %10.4f " % angles_s[0,e,k]
    stri += "\n"
fh.write(stri)
fh.close()
