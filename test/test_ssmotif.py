
import barnaba as bb

fname = "data/GNRA.pdb"
fname1 = "data/1S72.pdb"
dist = bb.ss_motif(fname,fname1,bulges=1,treshold=0.6,out='ss_motif')
fh=open("ssmotif.test.dat",'w')
stri = ""
for el in dist:
    stri += " %8.5f " % el[1]
    for jj in el[2]:
        stri += "%s-" % jj
    stri += "\n"
fh.write(stri)
fh.close()
