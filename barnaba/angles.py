import mdtraj as md
import numpy as np
import sys
import nucleic
import definitions

def backbone_angles_traj(traj):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    idxs = (nn.get_bb_torsion_idx()).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    angles = md.compute_dihedrals(traj,idxs,opt=False)
    # set to NaN where atoms are missing
    angles[:,missing[0]] = np.nan
    angles = angles.reshape((traj.n_frames,len(rna_seq),7))
    
    return rna_seq, angles
    
    
def backbone_angles(filename,topology=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    rna_seq,angles  = backbone_angles_traj(traj)
    return rna_seq,angles



def sugar_angles_traj(traj):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    idxs = (nn.get_sugar_torsion_idx()).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    angles = md.compute_dihedrals(traj,idxs,opt=False)
    # set to NaN where atoms are missing
    angles[:,missing[0]] = np.nan
    angles = angles.reshape((traj.n_frames,len(rna_seq),5))

    # sugar puckers
    x1 = angles[:,:,4] +  angles[:,:,1] -  angles[:,:,3] -   angles[:,:,0]
    x2 = 3.0776835*angles[:,:,2]
    phase = np.arctan2(x1,x2)
    phase[np.where(phase<0.0)] += 2.0*np.pi
    tm = angles[:,:,2]/np.cos(phase)
    angles2 = np.dstack((angles,phase,tm))

    return rna_seq, angles2
    
    
def sugar_angles(filename,topology=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    rna_seq,angles  = sugar_angles_traj(traj)
    return rna_seq,angles
