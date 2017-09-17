
import barnaba.ds_motif as bb

fname = "data/SARCIN.pdb"
fname1 = "data/1S72.pdb"

dist = bb.dsmotif(fname,fname1,l1=8,l2=7,bulges=0,treshold=0.65,write="test")
fh=open("dsmotif.dat",'w')
stri = ""
for el in dist:
    stri += " %8.5f " % el[1]
    for jj in el[2]:
        stri += "%s-" % jj
    stri += "\n"
fh.write(stri)
fh.close()
#fh = open("ermsd.dat",'w')
#stri = ""
#for e in range(len(dist)):
#    stri += "%14e \n" % (dist[e])

#fh.write(stri)
#fh.close()
