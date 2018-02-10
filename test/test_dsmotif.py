from __future__ import absolute_import, division, print_function
import barnaba as bb
import os
import glob
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))
fname = "%s/test/data/SARCIN.pdb" % cwd
fname1 = "%s/test/data/1S72.pdb" % cwd

def test_ssmotif():
    
    dist = bb.ds_motif(fname,fname1,l1=8,l2=7,bulges=0,threshold=0.65,out='%s/ds_motif' % outdir)
    stri = ""
    for el in dist:
        stri += " %8.5f " % el[1]
        for jj in el[2]:
            stri += "%s-" % jj
        stri += "\n"
        
    fh = open("%s/dsmotif_01.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/dsmotif_01.test.dat" % refdir)

    # compare all pdbs
    of = glob.glob('%s/ds_motif*.pdb' % refdir)
    for f in of:
        comp(f)

