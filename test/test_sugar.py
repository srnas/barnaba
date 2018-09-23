from __future__ import absolute_import, division, print_function
import barnaba as bb
import barnaba.definitions as dd
import os
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd
fname1 = "%s/test/data/samples.xtc" % cwd

def test_sugar_1():
    
    sugar_b,rr = bb.sugar_angles(fname)

    fh = open("%s/sugar_01.test.dat" % outdir,'w')
    stri = "%20s " % "#"
    for pp in dd.sugar_angles:
        stri += " %10s " % pp
    stri += "\n"
    
    for e in range(sugar_b.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(sugar_b.shape[2]):
            stri += " %10.4f " % sugar_b[0,e,k]
        stri += "\n"
    fh.write(stri)
    fh.close()
    comp("%s/sugar_01.test.dat" % refdir)

    

def test_sugar_2():
    
    resi = ["RC5_1_0","RG_69_0"]
    angles = ["nu1","nu4","nu3"]
    sugar_b,rr = bb.sugar_angles(fname,residues=resi,angles=angles)
    stri = "%20s " % "#"
    for pp in angles:
        stri += " %10s " % pp
    stri += "\n"

    for e in range(sugar_b.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(sugar_b.shape[2]):
            stri += " %10.4f " % sugar_b[0,e,k]
        stri += "\n"
        
    fh = open("%s/sugar_02.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/sugar_02.test.dat" % refdir)



def test_sugar_3():
    
    resi = ["RG_69_0","RU_37_0"]
    na = ["Phase","tm"]
    angles,rr = bb.pucker_angles(fname,residues=resi,altona=True)

    stri = "".join([" %10s " % pp for pp in na])
    stri += "\n"
    
    for e in range(angles.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(angles.shape[2]):
            stri += " %10.4f " % angles[0,e,k]
        stri += "\n"
    
    fh = open("%s/sugar_03.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/sugar_03.test.dat" % refdir)


def test_sugar_4():

    resi = ["RC5_1_0","RG_69_0","RU_37_0"]
    na = ["Phase","tm"]

    angles_b,rr = bb.pucker_angles(fname,residues=resi)
    stri = "%20s " % "#"
    for pp in na:
        stri += " %10s " % pp
    stri += "\n"

    for e in range(angles_b.shape[1]):
        stri += "%20s " % (rr[e])
        for k in range(angles_b.shape[2]):
            stri += " %10.4f " % angles_b[0,e,k]
        stri += "\n"

    fh = open("%s/sugar_04.test.dat" % outdir,'w')
    fh.write(stri)
    fh.close()
    comp("%s/sugar_04.test.dat" % refdir)
        

