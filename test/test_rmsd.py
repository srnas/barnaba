import barnaba as bb

# align pdb to pdb with the same sequence
fname = "data/sample1.pdb"
fname1 = "data/sample2.pdb"

dist = bb.rmsd(fname,fname1,out='aligned_1.pdb')
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("rmsd_01.test.dat",'w')
fh.write(stri)
fh.close()

# align pdb to pdb with different sequence
fname = "data/4v7t-pdb-bundle3_G521_00006.align.pdb"
fname1 = "data/centroid_10.pdb"

dist = bb.rmsd(fname,fname1,out='aligned_2.pdb')
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("rmsd_02.test.dat",'w')
fh.write(stri)
fh.close()

# align trajectory to pdb
fname = "data/sample1.pdb"
fname1 = "data/samples.xtc"

dist = bb.rmsd(fname,fname1,topology=fname,out='aligned_3.xtc')
stri = "".join([ "%14e \n" % (dd) for dd in dist])
fh = open("rmsd_03.test.dat",'w')
fh.write(stri)
fh.close()
