from scipy.spatial import distance
import definitions
import nucleic
import numpy as np
import mdtraj as md

###########################  functions #################

def calc_lcs(coords):
    # calculate center of mass
    origo = np.sum(coords,axis=0)/3.0
    
    # CoM-C2 (x axis)
    x = coords[0]-origo
    x_norm = np.sqrt(np.sum(x*x,axis=1))
    x = x/x_norm[:,np.newaxis]
    # CoM-C4/C6 
    c = coords[1]-origo
    # z/y axis
    z = np.cross(x,c,axis=1)
    z_norm = np.sqrt(np.sum(z*z,axis=1))
    z = z/z_norm[:,np.newaxis]
    
    y = np.cross(z,x,axis=1)
    lcs = np.array([x.T,y.T,z.T]).T
    return lcs, origo


def calc_3dmat(coords,cutoff):
    
    # prune search first
    lcs,origo= calc_lcs(coords)
    max_r  = np.max(definitions.f_factors)*cutoff
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.01))).T
        
    # calculate scaled distances
    #diff = [origo[y]-origo[x] for x,y in m_idx]
    diff = origo[m_idx[:,1]]-origo[m_idx[:,0]]
    
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])    
    return dotp,m_idx


def calc_rmat(coords,cutoff):

    dotp,m_idx = calc_3dmat(coords,cutoff)
    ll = coords.shape[1]
    
    mat = np.zeros((ll,ll,3))        
    
    if(dotp.shape[0]==0): return mat
    
    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    dotp[dotp_norm>cutoff] = 0.0
    
    
    mat[m_idx[:,0],m_idx[:,1]] = dotp
    # the matrix is not rescaled!
    return mat


def calc_mat_annotation(coords):
    
    lcs,origo= calc_lcs(coords)

    cutoff=1.58  # hardcoded 
    # prune search first
    max_r  = np.max(definitions.f_factors)*cutoff
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.001))).T
    
    if(len(m_idx)==0):
        return [],[]
    
    # calculate scaled distances
    #diff = [origo[y]-origo[x] for x,y in m_idx]
    diff = origo[m_idx[:,1]]-origo[m_idx[:,0]]
        
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_scale_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    angle = np.array([np.dot(lcs[i][:,2],lcs[j][:,2]) for i,j in  m_idx])

    ll = coords.shape[1]
        
    # create ll by ll matrix where ellipsoidal distance is less than cutoff
    cutoff_mat = np.zeros((ll,ll))
    cutoff_mat[m_idx[:,0],m_idx[:,1]] = dotp_scale_norm<cutoff
    
    #cutoff_mat[dotp_scale_norm<cutoff] = 1
    # symmetrize
    cutoff_mat *=cutoff_mat.T
    mat = np.zeros((ll,ll,3))
    mat[m_idx[:,0],m_idx[:,1]] = dotp
    
    mat *= cutoff_mat[:,:,np.newaxis]
    angles = np.zeros((ll,ll))
    angles[m_idx[:,0],m_idx[:,1]] = angle
    angles *= cutoff_mat
    return mat,angles

def calc_gmat(coords,cutoff):

 
    ll = coords.shape[1]

    mat = np.zeros((ll,ll,4))
    
    dotp,m_idx = calc_3dmat(coords,cutoff)
        
    # return zero matrix when there are no contacts
    if(dotp.shape[0]==0): return mat
    
    dotp *= np.array(definitions.scale)[np.newaxis,:]
    dotp_norm = np.sqrt(np.sum(dotp**2,axis=1))
    
    # calculate 4D g-vector
    ff = (np.pi*dotp_norm)/cutoff
    factor13 = np.sin(ff)/ff
    factor4= ((1.0+np.cos(ff))*cutoff)/np.pi
    gmat = dotp*factor13[:,np.newaxis]
    gmat = np.concatenate((gmat,factor4[:,np.newaxis]),axis=1)
    
    # set to zero when norm is larger than cutoff
    gmat[dotp_norm>cutoff] = 0.0
            
    mat[m_idx[:,0],m_idx[:,1]] = gmat
    
    return mat

###########################################
###########  ERMSD ########################

def ermsd_traj(reference,traj,cutoff=2.4):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = reference.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))
    
    coords_ref = reference.xyz[0,nn_ref.indeces_lcs]
    ref_mat = calc_gmat(coords_ref,cutoff).reshape(-1)
    #rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    gmats = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn_traj.indeces_lcs]
        gmats.append(calc_gmat(coords_lcs,cutoff).reshape(-1))
    dd = distance.cdist([ref_mat],gmats)/np.sqrt(len(nn_traj.ok_residues))
    return dd[0]

        
def dump_rvec_traj(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    rvecs = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        rvecs.append(calc_rmat(coords_lcs,cutoff))
        
    return np.asarray(rvecs), nn.ok_residues


def dump_gvec_traj(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    gvecs = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        gvecs.append(calc_gmat(coords_lcs,cutoff))
    return np.asarray(gvecs), nn.ok_residues


def backbone_angles_traj(traj,residues=None,angles=None):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =  nn.get_bb_torsion_idx(residues)
               
    if(angles==None):
        idx_angles = np.arange(all_idx.shape[1])
    else:
        idx_angles = [i for i in np.arange(all_idx.shape[1]) if(definitions.bb_angles[i]) in angles]
        assert(len(idx_angles)!=0)

    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    
    torsions = md.compute_dihedrals(traj,idxs,opt=False)
    
    # set to NaN where atoms are missing
    torsions[:,np.where(np.sum(idxs,axis=1)==0)[0]] = np.nan
    
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))
    
    return torsions, rr


def sugar_angles_traj(traj,residue=None):
    
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

    return angles2,nn.ok_residues



def rmsd_traj(reference,traj,out=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = reference.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)
    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))

    # loop over residues and find common heavy atoms
    idx_ref = []
    idx_target = []
    
    for ii in range(len(nn_ref.ok_residues)):
        
        res1 = nn_ref.ok_residues[ii]
        res2 = nn_traj.ok_residues[ii]

        resname1 = nn_ref.rna_seq_id[ii]
        resname2 = nn_traj.rna_seq_id[ii]
        
        # if the nucleotide is the same, use all atoms
        if(res1.name == res2.name):
            
            name2 = [at.name for at in res2.atoms if  at.name in definitions.nt_atoms[resname2]]
            #name1 = [at.name for at in res1.atoms if  at.name in definitions.nb_atoms[res1.name]]
            #print name2,name1
            for at in res1.atoms:
                if at.name in definitions.nt_atoms[resname1]:
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
        
    traj.superpose(reference,atom_indices=idx_target, ref_atom_indices=idx_ref)
    if(out!=None):
        traj.save(out)
    rmsd = np.sqrt(3*np.mean((traj.xyz[:, idx_target, :] - reference.xyz[0,idx_ref, :])**2, axis=(1,2)))
    return rmsd
