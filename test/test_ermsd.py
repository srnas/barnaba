import barnaba as bb
import os
import difflib
import sys

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))

def comp(s1,s2,j):


    diff = [x for x in difflib.unified_diff(s1,s2)]
    
    if(len(list(diff))==0):
        return 0
    else:
        #print '\n'.join(diff)
        sys.stderr.write(''.join(diff))
        fh = open("%s/diff_%d.test.dat" % (outdir,j),'w')
        fh.write(''.join(diff))
        fh.close()
        assert(1==2)
    
    
    
def test_ermsd_1():

    # align pdb to pdb with the same sequence
    fname = "%s/test/data/sample1.pdb" % cwd
    fname1 = "%s/test/data/sample2.pdb" % cwd
    
    dist = bb.ermsd(fname,fname1)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    
    fh = open("%s/ermsd_01.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()

    ref_string = (open("%s/ermsd_01.test.dat" % refdir)).readlines()
    tmp_string = (open("%s/ermsd_01.test.dat" % outdir)).readlines()

    return comp(ref_string,tmp_string,1)

    
def test_ermsd_2():
    # align pdb to pdb with different sequence
    fname = "%s/test/data/4v7t-pdb-bundle3_G521_00006.align.pdb" % cwd
    fname1 = "%s/test/data/centroid_10.pdb" % cwd
    
    dist = bb.ermsd(fname,fname1)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/ermsd_02.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    
    ref_string = (open("%s/ermsd_02.test.dat" % refdir)).readlines()
    tmp_string = (open("%s/ermsd_02.test.dat" % outdir)).readlines()

    return comp(ref_string,tmp_string,2)


def test_ermsd_3():
    
    fname = "%s/test/data/4v7t-pdb-bundle3_G521_00006.align.pdb" % cwd
    fname1 = "%s/test/data/centroid_10.pdb" % cwd
    
    dist = bb.ermsd(fname,fname1,cutoff=5.0)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/ermsd_03.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    
    ref_string = (open("%s/ermsd_03.test.dat" % refdir)).readlines()
    tmp_string = (open("%s/ermsd_03.test.dat" % outdir)).readlines()

    return comp(ref_string,tmp_string,3)


def test_ermsd_4():
    
    # align trajectory to pdb
    fname = "%s/test/data/sample1.pdb" % cwd
    fname1 = "%s/test/data/samples.xtc" % cwd
    
    dist = bb.ermsd(fname,fname1,topology=fname)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/ermsd_04.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()

    ref_string = (open("%s/ermsd_04.test.dat" % refdir)).readlines()
    tmp_string = (open("%s/ermsd_04.test.dat" % outdir)).readlines()

    return comp(ref_string,tmp_string,4)


