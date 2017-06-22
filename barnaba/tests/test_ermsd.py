
import barnaba.ermsd as bb

fname = "data/sample1.pdb"
fname1 = "data/samples.xtc"

dist = bb.ermsd(fname,fname1,topology=fname)


fh = open("ermsd.dat",'w')
stri = ""
for e in range(len(dist)):
    stri += "%14e \n" % (dist[e])

fh.write(stri)
fh.close()
