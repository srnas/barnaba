import barnaba.barnaba as bb

# align pdb to pdb with the same sequence
fname = "data/sample1.pdb"
fname1 = "data/sample2.pdb"

dist = bb.ermsd(fname,fname1)
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("ermsd_01.dat",'w')
fh.write(stri)
fh.close()

# align pdb to pdb with different sequence
fname = "data/4v7t-pdb-bundle3_G521_00006.align.pdb"
fname1 = "data/centroid_10.pdb"

dist = bb.ermsd(fname,fname1)
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("ermsd_02.dat",'w')
fh.write(stri)
fh.close()

fname = "data/4v7t-pdb-bundle3_G521_00006.align.pdb"
fname1 = "data/centroid_10.pdb"

dist = bb.ermsd(fname,fname1,cutoff=5.0)
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("ermsd_03.dat",'w')
fh.write(stri)
fh.close()


# align trajectory to pdb
fname = "data/sample1.pdb"
fname1 = "data/samples.xtc"

dist = bb.ermsd(fname,fname1,topology=fname)
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("ermsd_04.dat",'w')
fh.write(stri)
fh.close()
