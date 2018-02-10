from __future__ import absolute_import, division, print_function

import barnaba as bb
import os
import sys
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

    
def test_ermsd_1():

    # align pdb to pdb with the same sequence
    fname = "%s/test/data/sample1.pdb" % cwd
    fname1 = "%s/test/data/sample2.pdb" % cwd
    
    dist = bb.ermsd(fname,fname1)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    
    fh = open("%s/ermsd_01.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()

    comp("%s/ermsd_01.test.dat" % refdir)
    return 0

    
def test_ermsd_2():
    
    # align pdb to pdb with different sequence
    fname = "%s/test/data/4v7t-pdb-bundle3_G521_00006.align.pdb" % cwd
    fname1 = "%s/test/data/centroid_10.pdb" % cwd
    
    dist = bb.ermsd(fname,fname1)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/ermsd_02.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()

    comp("%s/ermsd_02.test.dat" % refdir)
    return 0

    


def test_ermsd_3():
    
    fname = "%s/test/data/4v7t-pdb-bundle3_G521_00006.align.pdb" % cwd
    fname1 = "%s/test/data/centroid_10.pdb" % cwd
    
    dist = bb.ermsd(fname,fname1,cutoff=5.0)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/ermsd_03.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    
    comp("%s/ermsd_03.test.dat" % refdir)
    return 0


def test_ermsd_4():
    
    # align trajectory to pdb
    fname = "%s/test/data/sample1.pdb" % cwd
    fname1 = "%s/test/data/samples.xtc" % cwd
    
    dist = bb.ermsd(fname,fname1,topology=fname)
    stri = "".join([ "%14e \n" % (dd) for dd in dist])
    fh = open("%s/ermsd_04.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()

    comp("%s/ermsd_04.test.dat" % refdir)
    return 0



test_ermsd_1()
