
import barnaba as bb

fname = "data/SARCIN.pdb"
fname1 = "data/1S72.pdb"

dist = bb.ds_motif(fname,fname1,l1=8,l2=7,bulges=0,treshold=0.65,out="ds_motif")
fh=open("dsmotif.test.dat",'w')
stri = ""
for el in dist:
    stri += " %8.5f " % el[1]
    for jj in el[2]:
        stri += "%s-" % jj
    stri += "\n"
fh.write(stri)
fh.close()
