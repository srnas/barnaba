import barnaba as bb
import mdtraj as md
import barnaba.smm as smm
import os
import filecmp

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir %s" % (outdir))


def test_smm_1():
    
    fname = "%s/test/data/UUCG.pdb" % cwd
    traj = "%s/test/data/UUCG.xtc" % cwd
    
    gvec,seq = bb.dump_gvec(traj,topology=fname)
    lent = gvec.shape[0]
    gvec = gvec.reshape(lent,-1)[::5]
    s = smm.SMM(gvec,eps=0.5)
    
