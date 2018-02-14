from __future__ import absolute_import, division, print_function
import barnaba.escore as escore
import os
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd
traj = "%s/test/data/samples.xtc" % cwd
ffname = "%s/test/data/1S72.pdb" % cwd

def test_score():
    
    # set "force-field"
    ee = escore.Escore([ffname])
    
    ss =  ee.score(traj,topology=fname)
    fh = open("%s/score_01.test.dat" % outdir,'w')    
    stri = " ".join(["%10.4e \n" % (ss[e]) for e in range(len(ss))])
    fh.write(stri)
    fh.close()
    comp("%s/score_01.test.dat" % refdir)

        

    
