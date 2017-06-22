
import barnaba.couplings as bb

fname = "data/sample1.pdb"

rr, angles_s = bb.couplings(fname)

fh = open("couplings.dat",'w')
stri = ""
for e in range(len(rr)):
    stri += "%20s " % (rr[e])
    for k in range(10):
        stri += " %10.4f " % angles_s[0,e,k]
    stri += "\n"
fh.write(stri)
fh.close()
