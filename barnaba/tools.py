import definitions
from scipy.spatial import distance
import numpy as np


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
