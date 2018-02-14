from __future__ import absolute_import, division, print_function
import barnaba as bb
import os
from comp_mine import comp
import glob

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))
fname = "%s/test/data/GNRA.pdb" % cwd
fname1 = "%s/test/data/1S72.pdb" % cwd

def test_ssmotif():
    
    dist = bb.ss_motif(fname,fname1,bulges=1,threshold=0.6,out='%s/ss_motif' % outdir)
    stri = ""
    for el in dist:
        stri += " %8.5f " % el[1]
        for jj in el[2]:
            stri += "%s-" % jj
        stri += "\n"
        
    fh = open("%s/ssmotif_01.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/ssmotif_01.test.dat" % refdir)
    
    # compare all pdbs
    of = glob.glob('%s/ss_motif*.pdb' % refdir)
    for f in of:
        comp(f)
        
