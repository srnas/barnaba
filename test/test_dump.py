import barnaba as bb
import os
import difflib
import sys

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd

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

    ref_string = (open("%s/dump_01.test.dat" % refdir)).readlines()
    tmp_string = (open("%s/dump_01.test.dat" % outdir)).readlines()
    comp(ref_string,tmp_string,1)
    
    ref_string = (open("%s/dump_02.test.dat" % refdir)).readlines()
    tmp_string = (open("%s/dump_02.test.dat" % outdir)).readlines()
    comp(ref_string,tmp_string,2)

    
