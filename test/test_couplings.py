from __future__ import absolute_import, division, print_function
import barnaba as bb
import os
import filecmp
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

top = "%s/test/data/sample1.pdb" % cwd
traj = "%s/test/data/samples.xtc" % cwd

def test_couplings_1():
    
    angles,rr = bb.jcouplings(traj,topology=top)

    fh = open("%s/couplings_01.test.dat" % outdir,'w')    
    for e in range(angles.shape[0]):
        stri = ""
        for k in range(angles.shape[1]):
            for l in range(angles.shape[2]):
                stri += " %10.4f " % angles[e,k,l]
            stri += "\n"
    stri += "\n"
    fh.write(stri)
    fh.close()
    comp("%s/couplings_01.test.dat" % refdir)


def test_couplings_2():
    
    angles,rr = bb.jcouplings(traj,topology=top,raw=True)
    fh = open("%s/couplings_02.test.dat" % outdir,'w')    
    for e in range(angles.shape[0]):
        stri = ""
        for k in range(angles.shape[1]):
            for l in range(angles.shape[2]):
                stri += " %10.4f " % angles[e,k,l]
            stri += "\n"
    fh.write(stri)
    fh.close()
    comp("%s/couplings_02.test.dat" % refdir)
