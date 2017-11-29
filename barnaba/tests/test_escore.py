
import barnaba.escore as escore


fname = "data/1S72.pdb"
#fname = "data/sample1.pdb"
ee = escore.Escore([fname])

fname = "data/samples.xtc"
tname = "data/sample1.pdb"
ss =  ee.score(fname,tname)
for i in ss:
    print i
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
