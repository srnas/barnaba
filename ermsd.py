import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist
import sys
import nucleic
import definitions

def ermsd_traj(ref,traj,cutoff=2.4):
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = ref.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))
    
    coords_ref = ref.xyz[0,nn_ref.indeces_lcs]
    ref_mat = nn_ref.get_gmat(coords_ref,cutoff).reshape(-1)
    #rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    gmats = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn_traj.indeces_lcs]
        gmats.append(nn_ref.get_gmat(coords_lcs,cutoff).reshape(-1))
    dd = cdist([ref_mat],gmats)/np.sqrt(len(nn_traj.ok_residues))
    return dd[0]
        

def ermsd(reference,target,cutoff=2.4,topology=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return ermsd_traj(ref,traj)
