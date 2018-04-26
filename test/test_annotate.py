from __future__ import absolute_import, division, print_function
import barnaba as bb
import os
import filecmp
from comp_mine import comp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

fname = "%s/test/data/sample1.pdb" % cwd
fname1 = "%s/test/data/samples.xtc" % cwd


def test_annotate_1():
    
    stackings, pairings, res = bb.annotate(fname1,topology=fname)
    dotbr,ss = bb.dot_bracket(pairings,res)

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



    
    fh = open("%s/dotbracket_01.test.dat" % outdir,'w')        
    stri = ""
    for k in range(len(pairings)):
        stri += " %06d " % k
        stri += "%s \n" % dotbr[k]
    fh.write(stri)
    fh.close()
    
    comp("%s/stackings_01.test.dat" % refdir)
    comp("%s/pairings_01.test.dat" % refdir)
    comp("%s/dotbracket_01.test.dat" % refdir)



def test_annotate_2():
    
    fname = "%s/test/data/1S72.pdb" % cwd
    stackings, pairings, res = bb.annotate(fname)
    dotbr,ss = bb.dot_bracket(pairings,res)
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

    fh = open("%s/dotbracket_02.test.dat" % outdir,'w')        
    fh.write(str(dotbr[0]))
    fh.close()

    comp("%s/stackings_02.test.dat" % refdir)
    comp("%s/pairings_02.test.dat" % refdir)
    comp("%s/dotbracket_02.test.dat" % refdir)
    


