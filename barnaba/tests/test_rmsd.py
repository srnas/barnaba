
import barnaba.rmsd as bb

fname = "data/sample1.pdb"
fname1 = "data/samples.xtc"

dist = bb.rmsd(fname,fname1,topology=fname,out='aligned.xtc')


fh = open("rmsd.dat",'w')
stri = ""
for e in range(len(dist)):
    stri += "%14e \n" % (dist[e])

fh.write(stri)
fh.close()
