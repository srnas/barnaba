from __future__ import absolute_import, division, print_function
import barnaba as bb
import os
from comp_mine import comp
import sys

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd


def test_dump():
    
    fh = open("%s/dump_01.test.dat" % outdir,'w')    
    rvecs,resi = bb.dump_rvec(fname,cutoff=1.7)        
    stri = ""
    for i1 in range(len(resi)):
        for i2 in range(len(resi)):
            if(sum((rvecs[0,i1,i2])**2) > 1.0e-5):
                stri += "%10s %10s %14e %14e %14e \n" % (resi[i1],resi[i2],rvecs[0,i1,i2,0],rvecs[0,i1,i2,1],rvecs[0,i1,i2,2])
    fh.write(stri)
    fh.close()


    gvecs,resi = bb.dump_gvec(fname)
    fh = open("%s/dump_02.test.dat" % outdir,'w')    
    stri = ""
    for i1 in range(len(resi)):
        for i2 in range(len(resi)):
            if(sum((gvecs[0,i1,i2])**2) > 1.0e-5):
                stri += "%10s %10s %14e %14e %14e %14e  \n" % (resi[i1],resi[i2],gvecs[0,i1,i2,0],gvecs[0,i1,i2,1],gvecs[0,i1,i2,2],gvecs[0,i1,i2,3])
    fh.write(stri)
    fh.close()

    comp("%s/dump_01.test.dat" % refdir)
    comp("%s/dump_02.test.dat" % refdir)
    

    
