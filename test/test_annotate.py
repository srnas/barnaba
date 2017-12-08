import barnaba as bb
import os
import filecmp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd
fname1 = "%s/test/data/samples.xtc" % cwd


def test_annotate_1():
    
    stackings, pairings, res = bb.annotate(fname1,topology=fname)
    fh = open("%s/stackings_01.test.dat" % outdir,'w')
    
    stri = "# STACKING \n"
    for k in range(len(stackings)):
        stri += "# frame %d \n" % k
        for e in range(len(stackings[k][0])):
            stri += "%15s " % (res[stackings[k][0][e][0]])
            stri += "%15s " % (res[stackings[k][0][e][1]])
            stri += " %4s \n" % (stackings[k][1][e])
        stri += "\n"
    fh.write(stri)
    fh.close()
    
    fh = open("%s/pairings_01.test.dat" % outdir,'w')
    stri = "# PAIRING \n"
    for k in range(len(pairings)):
        stri += "# frame %d \n" % k
        for e in range(len(pairings[k][0])):
            stri += "%15s " % (res[pairings[k][0][e][0]])
            stri += "%15s " % (res[pairings[k][0][e][1]])
            stri += " %4s \n" % (pairings[k][1][e])
    fh.write(stri)
    fh.close()

    assert(filecmp.cmp("%s/stackings_01.test.dat" % outdir,"%s/stackings_01.test.dat" % refdir)==True)
    assert(filecmp.cmp("%s/pairings_01.test.dat" % outdir,"%s/pairings_01.test.dat" % refdir)==True)


def test_annotate_2():
    
    fname = "%s/test/data/1S72.pdb" % cwd
    stackings, pairings, res = bb.annotate(fname)
    
    fh = open("%s/stackings_02.test.dat" % outdir,'w')
    stri = "# STACKING \n"
    for e in range(len(stackings[0][0])):
        stri += "%15s " % (res[stackings[0][0][e][0]])
        stri += "%15s " % (res[stackings[0][0][e][1]])
        stri += " %4s \n" % (stackings[0][1][e])
    fh.write(stri)
    fh.close()

    fh = open("%s/pairings_02.test.dat" % outdir,'w')
    stri = "# PAIRING \n"
    for e in range(len(pairings[0][0])):
        stri += "%15s " % (res[pairings[0][0][e][0]])
        stri += "%15s " % (res[pairings[0][0][e][1]])
        stri += " %4s \n" % (pairings[0][1][e])
    fh.write(stri)
    fh.close()

    assert(filecmp.cmp("%s/stackings_02.test.dat" % outdir,"%s/stackings_02.test.dat" % refdir)==True)
    assert(filecmp.cmp("%s/pairings_02.test.dat" % outdir,"%s/pairings_02.test.dat" % refdir)==True)



