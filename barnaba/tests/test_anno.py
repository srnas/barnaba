
import barnaba.annotate as bb

fname = "data/1S72.pdb"

aa,rr, anno = bb.annotate(fname)


fh = open("annotation.dat",'w')
stri = ""
for e in range(len(rr[0])):
    stri += "%15s " % (aa[rr[0][e][0]])
    stri += "%15s " % (aa[rr[0][e][1]])
    stri += " %4s \n" % anno[0][e]

fh.write(stri)
fh.close()
