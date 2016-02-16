import definitions
from scipy.spatial import distance
import numpy as np
import itertools as its
import re


def get_rna(topology):

    return [atom.index for atom in topology.atoms if ((atom.residue.name == 'U') or (atom.residue.name == 'rU') or (atom.residue.name == 'RU') or (atom.residue.name == 'RU5') or (atom.residue.name == 'RU3') or (atom.residue.name == 'U3') or (atom.residue.name == 'U5') or (atom.residue.name == 'C') or (atom.residue.name == 'rC') or (atom.residue.name == 'RC') or (atom.residue.name == 'RC5') or (atom.residue.name == 'RC3') or (atom.residue.name == 'C3') or (atom.residue.name == 'C5') or (atom.residue.name == 'G') or (atom.residue.name == 'rG') or (atom.residue.name == 'RG') or (atom.residue.name == 'RG5') or (atom.residue.name == 'RG3') or (atom.residue.name == 'G3') or (atom.residue.name == 'G5') or (atom.residue.name == 'A') or (atom.residue.name == 'rA') or (atom.residue.name == 'RA') or (atom.residue.name == 'RA5') or (atom.residue.name == 'RA3') or (atom.residue.name == 'A3') or (atom.residue.name == 'A5'))]

def get_lcs_idx(topology):

    rna_idx = get_rna(topology)
    C2_idx = [ii for ii in rna_idx if(topology.atom(ii).name=="C2")]
    C4_idx = [ii for ii in rna_idx if(topology.atom(ii).name=="C4")]
    C6_idx = [ii for ii in rna_idx if(topology.atom(ii).name=="C6")]
    lcs_idx = [C2_idx,C4_idx,C6_idx]
    
    assert(len(lcs_idx[0])== len(lcs_idx[1]))
    assert(len(lcs_idx[0])== len(lcs_idx[2]))
    for j in range(len(lcs_idx[0])):
        # sanity check 
        assert (topology.atom(lcs_idx[0][j]).residue == topology.atom(lcs_idx[1][j]).residue )
        assert (topology.atom(lcs_idx[0][j]).residue == topology.atom(lcs_idx[2][j]).residue )
        # swap! 
        if(topology.atom(lcs_idx[0][j]).residue.name in definitions.pur):
            lcs_idx[1][j],lcs_idx[2][j] = lcs_idx[2][j],lcs_idx[1][j]
            
    return lcs_idx

def get_other_idx(topology,atom):

    rna_idx = get_rna(topology)
    return [ii for ii in rna_idx if(topology.atom(ii).name==atom)]

def get_lcs(coords1,coords2,coords3):
    coords = np.array([coords1,coords2,coords3])
    
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
    return lcs,origo


def get_3dmat(coords1,coords2,coords3,cutoff):
    
    # get lcs
    lcs,origo = get_lcs(coords1,coords2,coords3)
    
    # prune search first
    max_r  = np.max(definitions.f_factors)*cutoff
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.01))).T
    
    # calculate scaled distances
    diff = [origo[y]-origo[x] for x,y in m_idx]
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
    return dotp,m_idx

def get_3dmat_square(coords1,coords2,coords3,cutoff):
    

    dotp,m_idx = get_3dmat(coords1,coords2,coords3,cutoff)
    if(dotp.shape[0]==0): return mat

    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    dotp[dotp_norm>cutoff] = 0.0
    ll = len(coords1)    
    mat = np.zeros((ll,ll,3))        
    mat[m_idx[:,0],m_idx[:,1]] = dotp
    # the matrix is not rescaled!
    return mat


def get_mat_annotation(coords1,coords2,coords3,cutoff):
    
    # get lcs
    lcs,origo = get_lcs(coords1,coords2,coords3)
    
    # prune search first
    max_r  = np.max(definitions.f_factors)*cutoff
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.001))).T

    if(len(m_idx)==0):
        return [],[]
    
    # calculate scaled distances
    diff = [origo[y]-origo[x] for x,y in m_idx]
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_scale_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    angle = np.array([np.dot(lcs[i][:,2],lcs[j][:,2]) for i,j in  m_idx])

    ll = len(coords1)
    
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
    
    return mat, angles

def get_gmat(coords1,coords2,coords3,cutoff):

    ll = len(coords1)
    mat = np.zeros((ll,ll,4))

    dotp,m_idx = get_3dmat(coords1,coords2,coords3,cutoff)
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
    
    #mat = np.zeros((ll,ll,4))
    mat[m_idx[:,0],m_idx[:,1]] = gmat
    
    return mat

def get_mat_score(coords1,coords2,coords3,cutoff):

    dotp,m_idx = get_3dmat(coords1,coords2,coords3,cutoff)
    if(dotp.shape[0]==0): return mat

    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    
    # set to zero when norm is larger than cutoff
    dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    dotp[dotp_norm>cutoff] = 0.0
    nonzero = dotp[~np.all(dotp == 0, axis=1)]
    
    return nonzero.T

def get_other_mat(coords1,coords2,coords3,coords4,cutoff):

    # get lcs
    lcs,origo = get_lcs(coords1,coords2,coords3)
    
    # prune search first
    max_r  = np.max(definitions.f_factors)*cutoff
    dmat = distance.cdist(origo,coords4)

    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.01))).T
    
    # calculate scaled distances
    diff = [coords4[y]-origo[x] for x,y in m_idx]
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    
    dotp[dotp_norm>cutoff] = 0.0
    
    mat = np.zeros((coords1.shape[0],coords4.shape[0],3))        
    mat[m_idx[:,0],m_idx[:,1]] = dotp
    return mat


def get_pattern(query):
    # build pattern for regular expression
    pattern = "^"
    for res in query:
        assert res in definitions.known_abbrev, "# Fatal error: character %s not known. Use AUCG/NYR" % (res)
        if(res in definitions.rna):
            pattern += res
        else:
            if(res == "N"):
                pattern += "[AUCG]"
            if(res == "Y"):
                pattern += "[UC]"
            if(res == "R"):
                pattern += "[AG]"
    pattern += "$"
    return pattern

def get_idx(sequence,query,bulges=0):

    ll = len(query)
    seq_str  = "".join(sequence)
    pattern = get_pattern(query)

    # generate first all possible indeces
    indeces = []
    for b in range(bulges+1):

        # create position of insertions
        comb= its.combinations(range(1,ll+b-1),b)
        for it1 in comb:
            idx1 = range(ll+b)
            # remove them from list
            for it2 in it1:
                idx1.remove(it2)
            for start in range(len(sequence)-ll-b+1):
                # get sequence
                idx2 =  [idx1[i] + start for i in range(ll)]
                substr = "".join([sequence[i] for i in idx2])
                if(re.match(pattern,substr) != None):
                    indeces.append(idx2)
    return indeces

def get_sup_idx(topology):

    rna_idx = get_rna(topology)
    return [ii for ii in rna_idx if(((definitions.pypu_dict[topology.atom(ii).residue.name] == "R") and topology.atom(ii).name in definitions.align_atoms_pur) or \
                                    ((definitions.pypu_dict[topology.atom(ii).residue.name] == "Y") and topology.atom(ii).name in definitions.align_atoms_pyr))]

        

def pdb2gmat(pdb,cutoff):
    
    idx = get_lcs_idx(pdb.topology)
    c1 = pdb.xyz[0,idx[0]]
    c2 = pdb.xyz[0,idx[1]]
    c3 = pdb.xyz[0,idx[2]]
    return get_gmat(c1,c2,c3,cutoff)


def pdb2gmatscore(pdb,cutoff):
    
    idx = get_lcs_idx(pdb.topology)
    c1 = pdb.xyz[0,idx[0]]
    c2 = pdb.xyz[0,idx[1]]
    c3 = pdb.xyz[0,idx[2]]
    return get_mat_score(c1,c2,c3,cutoff)



#def is_complementary(r1,r2):
#
#
#    if(defi)
#    if(r1=="A" and r2 == "U"): return True
#    if(r2=="A" and r1 == "U"): return True
#    if(r1=="G" and r2 == "C"): return True
#    if(r2=="G" and r1 == "C"): return True
#    return False
#

def wc_gaussian(vec):

     dev1  = vec - definitions.wc_mean
     maha1 = np.einsum('...k,...kl,...l->...', dev1, definitions.inv_sigma, dev1)
     return (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha1)
