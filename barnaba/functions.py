from scipy.spatial import distance
import sys
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

###########################################
###########  DUMP ########################
        
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


########### BACKBONE ANGLES ##############

def backbone_angles_traj(traj,residues=None,angles=None):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =  nn.get_bb_torsion_idx(residues)
               
    if(angles==None):
        idx_angles = np.arange(all_idx.shape[1])
    else:
  
        idx_angles = []
        for i in range(len(angles)):
            if(angles[i] in definitions.bb_angles):
                idx_angles.append(definitions.bb_angles.index(angles[i]))
            else:
                msg = "# Fatal error. requested angle \"%s\" not available.\n" % angles[i]
                msg += "# Choose from: %s \n" % definitions.bb_angles
                sys.stderr.write(msg)
                sys.exit(1)
   

    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    
    torsions = md.compute_dihedrals(traj,idxs,opt=False)
    
    # set to NaN where atoms are missing
    torsions[:,np.where(np.sum(idxs,axis=1)==0)[0]] = np.nan
    
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))
    
    return torsions, rr


############## couplings ##############

def jcouplings_traj(traj,residues=None,couplings=None,raw=False):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =nn.get_coupling_idx(residues)

    # if not provided, calculate all couplings
    if(couplings==None):
        idx_angles = np.arange(all_idx.shape[1])
        couplings = [str(el) for el in definitions.couplings_idx.keys()]
        idx_angles1 = [el for el in definitions.couplings_idx.values()]
    else:
        idx_angles = []
        idx_angles1 = []
        for i in range(len(couplings)):
            if(couplings[i] in definitions.couplings_idx.keys()):
                idx_angles.append(definitions.couplings_idx[couplings[i]])
                idx_angles1.append(i)
            else:
                msg = "# Fatal error. requested coupling \"%s\" not available.\n" % couplings[i]
                msg += "# Choose from: %s \n" % couplings_idx
                sys.stderr.write(msg)
                sys.exit(1)
   
    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    
    torsions = md.compute_dihedrals(traj,idxs,opt=False)
    
    # set to NaN where atoms are missing
    torsions[:,np.where(np.sum(idxs,axis=1)==0)[0]] = np.nan
    
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))

    #jcouplings = np.zeros((torsions.shape[0],torsions.shape[1],len(couplings)))
    # now calculate couplings
    if(raw):
        return torsions,rr

    jcouplings = np.empty((traj.n_frames,all_idx.shape[0],len(couplings)))*np.nan

    for i in range(len(couplings)):
        # get karplus coefficients
        
        coef = definitions.couplings_karplus[couplings[i]]
        #ii = definitions.couplings_idx[couplings[i]]
        ii =  idx_angles1[i]

        #print ii, couplings[i],
        angles = np.copy(torsions[:,:,ii])
        
        # add phase
        if(coef[4] != 0.0):
            angles += coef[4]
        cos = np.cos(angles)
        val = coef[0]*cos*cos + coef[1]*cos
        # add shift
        if(coef[2] != 0.0):
            val += coef[2]
        # add generalized karplus term
        if(coef[3] != 0.0):
            sin = np.sin(angles)
            val += coef[3]*cos*sin
        jcouplings[:,:,i] = val
    
    return jcouplings,rr

####################################################
############  SUGAR TORSION ########################

def sugar_angles_traj(traj,residues=None,angles=None):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =  nn.get_sugar_torsion_idx(residues)
    if(angles==None):
        idx_angles = np.arange(all_idx.shape[1])
    else:
        # find indeces corresponding to angles
        idx_angles = []
        for i in range(len(angles)):
            if(angles[i] in definitions.sugar_angles):
                idx_angles.append(definitions.sugar_angles.index(angles[i]))
            else:
                msg = "# Fatal error. requested angle \"%s\" not available.\n" % angles[i]
                msg += "# Choose from: %s \n" % (definitions.sugar_angles)
                sys.stderr.write(msg)
                sys.exit(1)

    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)

    torsions = md.compute_dihedrals(traj,idxs,opt=False)
    # set to NaN where atoms are missing
    torsions[:,missing[0]] = np.nan
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))

    return torsions, rr


########### PUCKER ANGLES ##############
def pucker_traj(traj,residues=None):

    torsions,rr = sugar_angles_traj(traj,residues=residues)
    x1 = torsions[:,:,4] +  torsions[:,:,1] -  torsions[:,:,3] -   torsions[:,:,0]
    x2 = 3.0776835*torsions[:,:,2]
    phase = np.arctan2(x1,x2)
    phase[np.where(phase<0.0)] += 2.0*np.pi
    tm = torsions[:,:,2]/np.cos(phase)
    angles = np.dstack((phase,tm))
    return angles, rr

######### scalar couplings #########


########### RMSD TRAJ  ##############
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



############### annotation  ##########
######## TODO ##############

def annotate_traj(traj):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)

    max_r  = np.max(definitions.f_factors)*1.58

    condensed_idx =  np.triu_indices(len(nn.ok_residues), 1)
    pairs = []
    annotations = []
    for i in xrange(traj.n_frames):

        # calculate LCS
        coords = traj.xyz[i,nn.indeces_lcs]
        lcs,origo= calc_lcs(coords)

        # normal distance
        dmat = distance.pdist(origo)
        # prune search
        m_idx = np.where(dmat<max_r)

        i1 = condensed_idx[0][m_idx]
        i2 = condensed_idx[1][m_idx]

        # calculate coordinates in LCS for short-ranged residues
        diff = origo[i1]-origo[i2]
        
        dotp_12 = np.asarray([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),i1)])
        dotp_21 = np.asarray([np.dot(-diff[i],lcs[j]) for i,j in zip(range(len(diff)),i2)])
        angle = np.asarray([np.dot(lcs[i][:,2],lcs[j][:,2]) for i,j in zip(i1,i2)])
        
        # zeta and rho 
        z_12 = dotp_12[:,2]**2
        z_21 = dotp_21[:,2]**2
        rho_12 = dotp_12[:,0]**2 + dotp_12[:,1]**2
        rho_21 = dotp_21[:,0]**2 + dotp_21[:,1]**2

        # find stacking here
        stackz_12 = np.where((z_12>0.04) & (z_12 <=0.225))
        stackz_21 = np.where((z_21>0.04) & (z_21 <=0.225))
        rhoz_12 =  np.where(rho_12<0.0625)
        rhoz_21 =  np.where(rho_21<0.0625)
        anglez = np.where(np.abs(angle) >0.76)
        #stackz_21, rhoz_12 , rhoz_21,anglez  
        # find pairing
        pairz_12 = np.where((z_12<=0.04))
        pairz_21 = np.where((z_21<=0.04))
        #print stackz_12
        #print stackz_21
        inter_zeta = np.intersect1d(stackz_12,stackz_21)
        union_rho = np.union1d(rhoz_12[0],rhoz_21[0])
        inter_stack = np.intersect1d(union_rho,inter_zeta)
        inter_stack = np.intersect1d(inter_stack,anglez)
        #print inter_stack
        #print anglez
        for p in inter_stack:
            print nn.rna_seq[i1[p]],
            print nn.rna_seq[i2[p]]
            
        #for ii in range(len(i1)):
        #    print nn.rna_seq[i1[ii]],
        #    print nn.rna_seq[i2[ii]],

        #    print dotp_12[ii],
        #    print dotp_21[ii],
        #    print np.sqrt(rho_12[ii]), np.sqrt(rho_21[ii]),
        #    print np.sqrt(z_12[ii]), np.sqrt(z_21[ii])
            #print i1[stackz_21],i2[stackz_21]
        exit()
        

############# cluster #############
