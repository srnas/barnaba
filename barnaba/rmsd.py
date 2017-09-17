import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist
import sys
import nucleic
import definitions

def rmsd_traj(ref,traj,out=None):
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = ref.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)
    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))

    # loop over residues and find common heavy atoms
    idx_ref = []
    idx_target = []
    for res1,res2 in zip(nn_ref.ok_residues,nn_traj.ok_residues):
        
        # if the nucleotide is the same, use all atoms
        if(res1.name == res2.name):
            name2 = [at.name for at in res2.atoms if  at.name in definitions.nt_atoms[res2.name]]
            #name1 = [at.name for at in res1.atoms if  at.name in definitions.nb_atoms[res1.name]]
            #print name2,name1
            for at in res1.atoms:
                if at.name in definitions.nt_atoms[res1.name]:
                    if(at.name in name2):
                        #print at.name
                        idx_ref.append(at.index)
                        idx_target.append(((res2.atom(at.name)).index))
        # else, use bb only
        else:
            name2 = [at.name for at in res2.atoms if  at.name in definitions.bb_atoms]
            for at in res1.atoms:
                if at.name in definitions.bb_atoms:
                    if(at.name in name2):
                        idx_ref.append(at.index)
                        idx_target.append(((res2.atom(at.name)).index))
            
    print "# found ",len(idx_ref), "atoms in common"
    
    if(len(idx_ref)<3):
        warn =  "# Only  %d atoms in common. abort.\n" % len(idx_ref)
        sys.stderr.write(warn)
        sys.exit(1)
        
    traj.superpose(ref,atom_indices=idx_target, ref_atom_indices=idx_ref)
    if(out!=None):
        traj.save(out)
    rmsd = np.sqrt(3*np.mean((traj.xyz[:, idx_target, :] - ref.xyz[0,idx_ref, :])**2, axis=(1,2)))
    return rmsd


def rmsd(reference,target,topology=None,out=None,quiet=True):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    if(quiet==False):
        sys.stderr.write(warn)

    return rmsd_traj(ref,traj,out)
