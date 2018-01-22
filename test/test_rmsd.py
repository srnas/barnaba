import barnaba as bb
import os
import filecmp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))

def test_rmsd_1():
    # align pdb to pdb with the same sequence
    fname = "%s/test/data/sample1.pdb" % cwd
    fname1 = "%s/test/data/sample2.pdb" % cwd

    dist = bb.rmsd(fname,fname1,out='%s/aligned_1.pdb' % outdir)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    
    fh = open("%s/rmsd_01.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    assert(filecmp.cmp("%s/rmsd_01.test.dat" % outdir,"%s/rmsd_01.test.dat" % refdir)==True)
    #assert(filecmp.cmp("%s/aligned_1.pdb" % outdir,"%s/aligned_1.pdb" % refdir)==True)

def test_rmsd_2():
    # align pdb to pdb with different sequence
    fname = "%s/test/data/4v7t-pdb-bundle3_G521_00006.align.pdb" % cwd
    fname1 = "%s/test/data/centroid_10.pdb" % cwd

    dist = bb.rmsd(fname,fname1,out='%s/aligned_2.pdb' % outdir)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/rmsd_02.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    assert(filecmp.cmp("%s/rmsd_02.test.dat" % outdir,"%s/rmsd_02.test.dat" % refdir)==True)
    #assert(filecmp.cmp("%s/aligned_2.pdb" % outdir,"%s/aligned_2.pdb" % refdir)==True)
    
def test_rmsd_3():
    
    # align trajectory to pdb
    fname = "%s/test/data/sample1.pdb" % cwd
    fname1 = "%s/test/data/samples.xtc" % cwd

    
    dist = bb.rmsd(fname,fname1,topology=fname,out='%s/aligned_3.xtc' % outdir)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])

    fh = open("%s/rmsd_03.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    assert(filecmp.cmp("%s/rmsd_03.test.dat" % outdir,"%s/rmsd_03.test.dat" % refdir)==True)
    #assert(filecmp.cmp("%s/aligned_3.xtc" % outdir,"%s/aligned_3.xtc" % refdir)==True)
