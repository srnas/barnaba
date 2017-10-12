import mdtraj as md
import numpy as np
import sys
import nucleic
import definitions


def couplings_traj(traj,raw_angles=False):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    # find indeces
    idxs = (nn.get_couplings_torsion_idx()).reshape(-1,4)

    rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    angles = md.compute_dihedrals(traj,idxs,opt=False)

    # set to NaN where atoms are missing
    missing = np.where(np.sum(idxs,axis=1)==0)
    angles[:,missing[0]] = np.nan
    
    angles = angles.reshape((traj.n_frames,len(rna_seq),len(definitions.j3)))
    if(raw_angles):
        return rna_seq, angles

    # calculate couplings
    par = np.array([x[2] for x in definitions.j3])
    dd = angles + par[:,4][np.newaxis,np.newaxis,:]

    cos = np.cos(dd)
    sin = np.sin(dd)
    #print cos*cos*par[:,0][np.newaxis,np.newaxis,:]

    couplings = cos*cos*(par[:,0][np.newaxis,np.newaxis,:]) + \
                cos*(par[:,1][np.newaxis,np.newaxis,:]) + \
                par[:,2][np.newaxis,np.newaxis,:] + \
                cos*sin*par[:,3][np.newaxis,np.newaxis,:]
    
    return rna_seq, couplings
    
    
def couplings(filename,topology=None,raw_angles=False):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    rna_seq,couplings  = couplings_traj(traj,raw_angles=raw_angles)
    return rna_seq,couplings



def which_couplings():
    labels = [s[0] for s in definitions.j3]
    return labels
    
    
