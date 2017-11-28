from scipy.spatial import distance
import sys
import definitions
import nucleic
import numpy as np
import mdtraj as md
import itertools 


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


    cutoff_sq=2.89  # hardcoded cutoff squared  (1.7)
    # prune search first
    max_r  = np.max(definitions.f_factors)*np.sqrt(cutoff_sq)
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.001))).T
    
    if(len(m_idx)==0):
        return [],[]
    
    # calculate scaled distances
    diff = origo[m_idx[:,1]]-origo[m_idx[:,0]]
        
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_scale_norm_square = np.sum(dotp_scale**2,axis=1)

    # find pairs with low ellipsoidal distance
    low_idx = np.where(dotp_scale_norm_square<cutoff_sq)
    pairs_tmp = m_idx[low_idx]
    pairs_labs = ["%d_%d" % (aa[0],aa[1]) for aa in pairs_tmp]
    # remove cases in which r_ij < cutoff, r_ja > cutoff
    pairs = []
    vectors = []
    angles = []
    for i,aa in enumerate(pairs_tmp):
        rev = "%d_%d" % (aa[1],aa[0])
        if( (rev in pairs_labs) and (aa[1]>aa[0]) ):
            pairs.append([aa[0],aa[1]])
            other_idx = pairs_labs.index(rev)
            vectors.append([dotp[low_idx[0][i]], dotp[low_idx[0][other_idx]]])
            angles.append(np.dot(lcs[aa[0]][:,2],lcs[aa[1]][:,2]))

    return np.array(pairs), np.array(vectors), np.array(angles)
    
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
        
    return np.asarray(rvecs), nn.rna_seq


def dump_gvec_traj(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    gvecs = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        gvecs.append(calc_gmat(coords_lcs,cutoff))
    return np.asarray(gvecs), nn.rna_seq


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
def pucker_angles_traj(traj,residues=None):

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


### Single stranded motifs ##
def ss_motif_traj(ref,traj,treshold=0.8,cutoff=2.4,sequence=None,bulges=0,out=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = ref.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    ll = len(nn_ref.ok_residues)
    if(sequence==None):
        sequence = "N"*ll
    else:
        assert(len(sequence)==ll)
        
    coords_ref = ref.xyz[0,nn_ref.indeces_lcs]
    ref_mat = calc_gmat(coords_ref,cutoff).reshape(-1)
    
    rna_seq = nn_traj.rna_seq_id
    res_idxs = definitions.get_idx(rna_seq,sequence,bulges)
    resname_idxs = [[nn_traj.rna_seq[l] for l in rr]  for rr  in res_idxs]
    if(len(res_idxs)==0):
        return []

    lcs_idx = nn_traj.indeces_lcs
    results = []
    count = 1
    for i in xrange(traj.n_frames):
        
        gmats = [ calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1)  for j in res_idxs]
            
        dd = distance.cdist([ref_mat],gmats)/np.sqrt(ll)
        low = np.where(dd[0]<treshold)
        for k in low[0]:
            results.append([i,dd[0,k],resname_idxs[k]])

            # Write aligned PDB 
            if(out != None):
                pdb_out = "%s_%05d_%s_%d.pdb" % (out,count,resname_idxs[k][0],i)
                # slice trajectory
                tmp_atoms = []
                tmp_res =[]
                for r1 in res_idxs[k]:
                    tmp_atoms.extend([at.index for at in nn_traj.ok_residues[r1].atoms])
                    tmp_res.append([at for at in nn_traj.ok_residues[r1].atoms])
                traj_slice = traj[i].atom_slice(tmp_atoms)
                
                # align whatever is in common in the backbone
                idx_target = []
                idx_ref = []
                for res1,res2 in zip(nn_ref.ok_residues,traj_slice.topology.residues):
                    name2 = [at.name for at in res2.atoms if  at.name in definitions.bb_atoms]
                    for at in res1.atoms:
                        if at.name in definitions.bb_atoms:
                            if(at.name in name2):
                                idx_ref.append(at.index)
                                idx_target.append(((res2.atom(at.name)).index))
                                #idx_target.append(res2[name2.index(at.name)].index)
                traj_slice.superpose(ref,atom_indices=idx_target, ref_atom_indices=idx_ref)
                traj_slice.save(pdb_out)

            count += 1
    
    return results



##############################################################
###### DOuble stranded motifs ############


def ds_motif_traj(ref,traj,l1,l2,treshold=0.9,cutoff=2.4,sequence=None,bulges=0,out=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = ref.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    ll = len(nn_ref.ok_residues)
    assert(ll==l1+l2)
    if(sequence==None):
        sequence1 = "N"*l1
        sequence2 = "N"*l2
    else:
        sequence1=(sequence).split("%")[0]
        sequence2=(sequence).split("%")[1]
        assert len(sequence1)==l1, "# FATAL: query structure and sequence length mismatch!"
        assert len(sequence2)==l2, "# FATAL: query structure and sequence length mismatch!"

        
    coords_ref = ref.xyz[0,nn_ref.indeces_lcs]
    ref_mat = calc_gmat(coords_ref,cutoff).reshape(-1)

    coords_ref1 = ref.xyz[0,nn_ref.indeces_lcs[:,0:l1]]

    ref_mat1 = calc_gmat(coords_ref1,cutoff).reshape(-1)
    coords_ref2 = ref.xyz[0,nn_ref.indeces_lcs[:,l1:]]
    ref_mat2 = calc_gmat(coords_ref2,cutoff).reshape(-1)
   
    # calculate center of mass distances
    # this will be used to prune the search!
    ref_com1 = np.average(np.average(coords_ref1,axis=0),axis=0)
    ref_com2 = np.average(np.average(coords_ref2,axis=0),axis=0)
    dcom =  np.sqrt(np.sum((ref_com1-ref_com2)**2))
    
    # find indeces of residues according to sequence
    rna_seq = nn_traj.rna_seq_id

    all_idx1 = definitions.get_idx(rna_seq,sequence1,bulges)
    if(len(all_idx1)==0): return []
    all_idx2 = definitions.get_idx(rna_seq,sequence2,bulges)
    if(len(all_idx2)==0): return []


    lcs_idx = nn_traj.indeces_lcs
    idxs_combo = []
    results = []
    count = 1
    for i in xrange(traj.n_frames):

        # calculate eRMSD for strand1 
        gmats1 = [calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1) for j in all_idx1]
        dd1 = distance.cdist([ref_mat1],gmats1)
        low1 = np.where(dd1[0]<treshold*np.sqrt(l1))
        
        # calculate eRMSD for strand2 
        gmats2 = [calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1) for j in all_idx2]
        dd2 = distance.cdist([ref_mat2],gmats2)
        low2 = np.where(dd2[0]<treshold*np.sqrt(l2))

        # do combination
        combo_list = list(itertools.product(low1[0],low2[0]))
        
        gmats_combo = []
        
        for cc in combo_list:
            llc = len(set(all_idx1[cc[0]] + all_idx2[cc[1]]))
            # skip overlapping
            if(llc != l1 + l2): continue
            
            # skip distant
            com1 = np.average(np.average(traj.xyz[i,lcs_idx[:,all_idx1[cc[0]]]],axis=0),axis=0)
            com2 = np.average(np.average(traj.xyz[i,lcs_idx[:,all_idx2[cc[1]]]],axis=0),axis=0)
            dcoms = np.sqrt(np.sum((com1-com2)**2))
            if(dcoms > 2.5*dcom): continue

            idx_combo = all_idx1[cc[0]] + all_idx2[cc[1]]
            idxs_combo.append(idx_combo)
            gmats_combo.append(calc_gmat(traj.xyz[i,lcs_idx[:,idx_combo]],cutoff).reshape(-1))

        # calculate distances
        dd_combo = distance.cdist([ref_mat],gmats_combo)
        low_combo = np.where(dd_combo[0]<treshold*np.sqrt(l1 + l2))

        for k in low_combo[0]:
            
            #print idxs_combo[k]
            resname_idxs = [nn_traj.rna_seq[l] for l  in idxs_combo[k]]
            
            results.append([i,dd_combo[0,k]/np.sqrt(l1 + l2),resname_idxs])

            #print results[-1]
            # Write aligned PDB 
            if(out != None):
                pdb_out = "%s_%05d_%s_%d.pdb" % (out,count,resname_idxs[k][0],i)
                # slice trajectory
                tmp_atoms = []
                tmp_res =[]
                for r1 in idxs_combo[k]:
                    tmp_atoms.extend([at.index for at in nn_traj.ok_residues[r1].atoms])
                    tmp_res.append([at for at in nn_traj.ok_residues[r1].atoms])
                traj_slice = traj[i].atom_slice(tmp_atoms)
                
                # align whatever is in common in the backbone
                idx_target = []
                idx_ref = []
                for res1,res2 in zip(nn_ref.ok_residues,traj_slice.topology.residues):
                    name2 = [at.name for at in res2.atoms if  at.name in definitions.bb_atoms]
                    for at in res1.atoms:
                        if at.name in definitions.bb_atoms:
                            if(at.name in name2):
                                idx_ref.append(at.index)
                                idx_target.append(((res2.atom(at.name)).index))
                                #idx_target.append(res2[name2.index(at.name)].index)
                traj_slice.superpose(ref,atom_indices=idx_target, ref_atom_indices=idx_ref)
                traj_slice.save(pdb_out)

            count += 1
    return results

############### annotation  ##########
def wc_gaussian(vec):

    vec = 10.0*vec
    dev1  = vec - definitions.wc_mean
    maha1 = np.einsum('...k,...kl,...l->...', dev1, definitions.inv_sigma, dev1)
    return (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha1)


def dihedral(p1,p2,p3,p4):

    # difference vectors b0 is reversed
    #b0 = vecs[:,1] - vecs[:,0]
    #b1 = vecs[:,1] - vecs[:,2]
    #b2 = vecs[:,2] - vecs[:,3]
    b0 = p2-p1
    b1 = p2-p3
    b2 = p3-p4
    # norm
    norm_sq = np.sum(b1**2)
    norm_sq_inv = 1.0/norm_sq

    #print (np.sum(b0*b1,axis=1)*b1).shape
    v0 = b0 - b1*((np.sum(b0*b1)*norm_sq_inv))
    v2 = b2 - b1*((np.sum(b0*b2)*norm_sq_inv))
    x = np.sum(v0*v2)
    m = np.cross(v0,b1)*np.sqrt(norm_sq_inv)
    y = np.sum(m*v2)
    return np.arctan2( y, x )


def annotate_traj(traj):

    # this is the binning for annotation
    bins = [0,1.84,3.84,2.*np.pi]
    bins_label = ["W","H","S"]
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    
    max_r  = np.max(definitions.f_factors)*1.58
    condensed_idx =  np.triu_indices(len(nn.ok_residues), 1)

    stackings = []
    pairings = []
    others = []
    
    for i in xrange(traj.n_frames):

        # calculate LCS
        coords = traj.xyz[i,nn.indeces_lcs]

        # find bases in close contact (within ellipsoid w radius sqrt(2.5))
        pairs,vectors,angles = calc_mat_annotation(coords)

        # calculate rho
        rho_12 = vectors[:,0,0]**2 + vectors[:,0,1]**2
        rho_21 = vectors[:,1,0]**2 + vectors[:,1,1]**2

        # calculate z squared 
        z_12 = vectors[:,0,2]**2
        z_21 = vectors[:,1,2]**2
        # find stacked bases 

        # z_ij  AND z_ji > 2 AA
        stackz_12 = np.where(z_12>0.04)
        stackz_21 = np.where(z_21>0.04)
        
        # rho_ij OR rho_ji < 2.5 AA
        rhoz_12 =  np.where(rho_12<0.0625)
        rhoz_21 =  np.where(rho_21<0.0625)
        union_rho = np.union1d(rhoz_12[0],rhoz_21[0])

        # angle between normal planes < 40 deg
        anglez = np.where(np.abs(angles) >0.766)

        # intersect all criteria
        inter_zeta = np.intersect1d(stackz_12,stackz_21)
        inter_stack = np.intersect1d(union_rho,inter_zeta)
        stacked = np.intersect1d(inter_stack,anglez)

        # now I have found all stackings.
        #stacked_pairs = [[nn.rna_seq[pairs[k,0]],nn.rna_seq[pairs[k,1]]] for k in stacked]
        stacked_pairs = [[pairs[k,0],pairs[k,1]] for k in stacked]
        stacked_annotation = np.chararray((len(stacked_pairs),2))
        stacked_annotation[:] = ">"
        
        # revert where z_ij is negative
        rev1 = np.where(vectors[stacked,0,2]<0)
        stacked_annotation[rev1[0],0] = "<"

        rev2 = np.where(vectors[stacked,1,2]>0)
        stacked_annotation[rev2[0],1] = "<"

        stacked_annotation = ["".join(el) for el in stacked_annotation]
        ######################################################

        # find paired bases  (z_ij < 2 AA AND z_j1 < 2 AA)
        #paired = [j for j in xrange(len(pairs)) if((j not in stackz_12[0]) and (j not in stackz_21[0]))]
        paired = [j for j in xrange(len(pairs)) if((j not in stackz_12[0]) or (j not in stackz_21[0]))]
        #paired_pairs = [[nn.rna_seq[pairs[k,0]],nn.rna_seq[pairs[k,1]]] for k in paired]
        paired_pairs = [[pairs[k,0],pairs[k,1]] for k in paired]
        
        # calculate edge angle. subtract 0.16 as Watson edge is not zero
        edge_angles_1 = np.arctan2(vectors[paired,0,1],vectors[paired,0,0]) - definitions.theta1
        edge_angles_2 = np.arctan2(vectors[paired,1,1],vectors[paired,1,0]) - definitions.theta1

        # shift to 0-2pi range
        edge_angles_1[np.where(edge_angles_1<0.0)] += 2.*np.pi
        edge_angles_2[np.where(edge_angles_2<0.0)] += 2.*np.pi

        # find edge: 0 Watson, 1:Hoogsteen, 2sugar
        edge_1 = np.digitize(edge_angles_1,bins)-1
        edge_2 = np.digitize(edge_angles_2,bins)-1

        # calculate dihedral angle for cis/trans
        #angle_glyco = np.array([dihedral(traj.xyz[i,nn.indeces_glyco[1,pairs[k,0]]],\
        #                        traj.xyz[i,nn.indeces_glyco[0,pairs[k,0]]],\
        #                        traj.xyz[i,nn.indeces_glyco[0,pairs[k,1]]],\
        #                        traj.xyz[i,nn.indeces_glyco[1,pairs[k,1]]]) for k in paired])
    
        
        paired_annotation = np.chararray((len(paired_pairs),3))
        paired_annotation[:] = "X"
        
        # now explicit loop, a bit messy. sorry, Guido.
        for j in xrange(len(paired_pairs)):

            # if angle is larger than 60 deg, skip
            if(np.abs(angles[paired[j]]) < 0.5 ): continue

            index1 = paired_pairs[j][0]
            index2 = paired_pairs[j][1]
            
            # find donor and acceptor in first base
            r1_donor = nn.donors[index1]
            r1_acceptor = nn.acceptors[index1]
            # find donor and acceptor in second base
            r2_donor = nn.donors[index2]
            r2_acceptor = nn.acceptors[index2]
            combo_list = list(itertools.product(r1_donor,r2_acceptor)) + list(itertools.product(r1_acceptor,r2_donor))
            
            # distances between donor and acceptors
            delta = np.diff(traj.xyz[i,combo_list],axis=1)
            dist_sq = np.sum(delta**2,axis=2)
            # number or distances less than 3.3 AA is n_hbonds
            n_hbonds = (dist_sq<0.1089).sum()
            # if no hydrogen bonds, skip 
            if(n_hbonds==0): continue

            # find edge
            paired_annotation[j][0] = bins_label[edge_1[j]]
            paired_annotation[j][1] = bins_label[edge_2[j]]

            gidxs = [nn.indeces_glyco[index1][1],nn.indeces_glyco[index1][0],nn.indeces_glyco[index2][0],nn.indeces_glyco[index2][1]]

            # if atoms in glyco are missing, do not calculate cis/trans
            if(None in gidxs):
                paired_annotation[j][2] = "x"
            else:
                angle_glyco = dihedral(traj.xyz[i,gidxs[0]],traj.xyz[i,gidxs[1]],traj.xyz[i,gidxs[2]],traj.xyz[i,gidxs[3]])
                if(np.abs(angle_glyco) > 0.5*np.pi):
                    paired_annotation[j][2] = "t"
                else:
                    paired_annotation[j][2] = "c"
                
            # if is WWc, check for Watson-crick and GU
            if("".join(paired_annotation[j]) == "WWc"):
                r1 =  nn.rna_seq_id[index1]
                r2 =  nn.rna_seq_id[index2]
                ll = "".join(sorted([r1,r2]))
                if((ll=="AU" and n_hbonds > 1) or ( ll == "CG" and n_hbonds > 2)):
                    paired_annotation[j] = ["W","C","c"]
                if(ll=="GU" and n_hbonds > 1):
                    paired_annotation[j] = ["G","U","c"]
                    
            #print nn.rna_seq[index1],
            #print nn.rna_seq[index2],
            #print "".join(paired_annotation[j]),
            #print vectors[paired[j],0],
            #print vectors[paired[j],1],
            #print  angles[paired[j]]
        #exit()
        paired_annotation = ["".join(el) for el in paired_annotation]

        #unassigned = [j for j in xrange(len(pairs)) if(j not in stacked and j not in paired)]
        #unassigned_pairs = [[pairs[k,0],pairs[k,1]] for k in unassigned]
        #for k in unassigned:
        #    print nn.rna_seq[pairs[k,0]],nn.rna_seq[pairs[k,1]], vectors[k,0], vectors[k,1]
        
        stackings.append([stacked_pairs,stacked_annotation])
        pairings.append([paired_pairs,paired_annotation])
        #others.append(unassigned_pairs)

    return stackings, pairings, nn.rna_seq

############# cluster #############
