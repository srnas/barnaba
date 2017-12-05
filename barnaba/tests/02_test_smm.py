import mdtraj as md
import barnaba.smm as smm
import barnaba.barnaba as bb

fname = "data/sample1.pdb"

# initialize class. Only PDB are accepted
traj = "data/UUCG.xtc"
top = "data/UUCG.pdb"
gvec,seq = bb.dump_gvec(traj,top)
lent = gvec.shape[0]
gvec = gvec.reshape(lent,-1)[::5]

smm = smm.SMM(gvec,eps=0.5)

