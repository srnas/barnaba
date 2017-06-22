
import barnaba.ss_motif as bb

fname = "data/GNRA.pdb"
fname1 = "data/1S72.pdb"

dist = bb.ssmotif(fname,fname1,bulges=1,treshold=0.6,write='prova')
fh=open("ssmotif.dat",'w')
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
