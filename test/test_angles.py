from __future__ import absolute_import, division, print_function
import barnaba as bb
import barnaba.definitions as dd
import os
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd
fname1 = "%s/test/data/samples.xtc" % cwd

def test_angles_1():
    
    angles_b,rr = bb.backbone_angles(fname)

    fh = open("%s/angles_01.test.dat" % outdir,'w')
    stri = "%20s " % "#"
    for pp in dd.bb_angles:
        stri += " %10s " % pp
    stri += "\n"
    
    for e in range(angles_b.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(angles_b.shape[2]):
            stri += " %10.4f " % angles_b[0,e,k]
        stri += "\n"
    fh.write(stri)
    fh.close()
    comp("%s/angles_01.test.dat" % refdir)
    

def test_angles_2():
    
    resi = ["RC5_1_0","RG_69_0"]
    angles = ["alpha","beta","chi"]
    angles_b,rr = bb.backbone_angles(fname,residues=resi,angles=angles)
    stri = "%20s " % "#"
    for pp in angles:
        stri += " %10s " % pp
    stri += "\n"

    for e in range(angles_b.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(angles_b.shape[2]):
            stri += " %10.4f " % angles_b[0,e,k]
        stri += "\n"
        
    fh = open("%s/angles_02.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()

    comp("%s/angles_02.test.dat" % refdir)


def test_angles_3():
    
    resi = ["RG_69_0","RU_37_0"]
    angles = ["gamma","alpha"]
    
    angles_b,rr = bb.backbone_angles(fname,residues=resi,angles=angles)
    stri = "%20s " % "#"
    for pp in angles:
        stri += " %10s " % pp
    stri += "\n"

    for e in range(angles_b.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(angles_b.shape[2]):
            stri += " %10.4f " % angles_b[0,e,k]
        stri += "\n"

    fh = open("%s/angles_03.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/angles_03.test.dat" % refdir)



def test_angles_4():
    
    resi = ["RG_69_0","RU_37_0"]
    angles = ["gamma","alpha"]
    
    angles_b,rr = bb.backbone_angles(fname1,topology=fname,residues=resi,angles=angles)
    stri = ""
    for p in range(angles_b.shape[0]):
        for k in range(angles_b.shape[2]):
            stri += " %10.4f %10.4f " % (angles_b[p,0,k],angles_b[p,1,k])
        stri += "\n"
        
    fh = open("%s/angles_04.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/angles_04.test.dat" % refdir)
        


