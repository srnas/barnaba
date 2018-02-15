#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2017 Sandro Bottaro (sandro.bottaro@bio.ku.dk)
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import, division, print_function
from scipy.spatial import distance
from . import definitions
import numpy as np


def calc_lcs(coords):    
    """
    Calculate local coordinates system

    Calculate origin position and xyz orthonormal vectors in the six-member rings

    Parameters
    ----------
    coords : (3,n,3) numpy array 
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    Returns
    -------
    lcs : (n,3,3) numpy array
        x y z vectors defining a local coordinate system in the geometrical center of the nucleobase
    origo : (n,3) numpy array 
        origin  of local coordinate systems        
    """
    
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
    """
    Calculate relative position of nucleobases within an ellipsoidal cutoff

    Parameters
    ----------
    coords : (3,n,3) numpy array 
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    cutoff : float
       ellipsoidal cutoff
    
    Returns
    -------
    dotp : (x,3) numpy array
       xyz coordinates for each pair

    m_idx : (x,2) numpy array 
       indeces of the pair
    """
    
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
    """
    Calculate relative position of nucleobases. Returns a square matrix with zero elements for pairs outside ellispoidal cutoff 

    Parameters
    ----------
    coords : (3,n,3) numpy array 
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    cutoff : float
        ellipsoidal cutoff
    
    Returns
    -------
    dotp : (n,n,3) numpy array
       xyz coordinates for each pair. For pairs outside the cutoff the coordinates are (0,0,0)
    """

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

    """
    Calculate matrix for annotation purposes.

    Parameters
    ----------
    coords : (3,n,3) numpy array 
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    Returns
    -------
    pairs : (x,2) numpy array
        indeces of pairs to be considered when performing the annotation

    vectors: (x,2,3) numpy array
        position vector r_ij and r_ji 
    
    angles: (x) numpy array
       cosine of angle beween the normal vectors constructed on base i and base j
    """

    lcs,origo= calc_lcs(coords)

    cutoff_sq=2.89  # hardcoded cutoff squared  (1.7)
    # prune search first
    max_r  = np.max(definitions.f_factors)*np.sqrt(cutoff_sq)
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.001))).T
    
    if(len(m_idx)==0):
        return [],[], []
    
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
    
    """
    Calculate G-vectors for each pair of bases within ellipsoidal cutoff distance 

    Parameters
    ----------
    coords : (3,n,3) numpy array 
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    cutoff : float
        ellipsoidal cutoff
    
    Returns
    -------
    dotp : (n,n,4) numpy array
        G coordinates for each pair. For pairs outside the cutoff the coordinates are (0,0,0,0)
    """

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

def calc_scoremat(coords,cutoff):
    """
    Calculate relative position of nucleobases within an ellipsoidal cutoff

    Parameters
    ----------
    coords : (3,n,3) numpy array 
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    cutoff : float
        ellipsoidal cutoff
    
    Returns
    -------
    dotp : (x,3) numpy array
        xyz coordinates for each pair
    """
    
    dotp,m_idx = calc_3dmat(coords,cutoff)
    ll = coords.shape[1]
        
    dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
    dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    dotp[dotp_norm>cutoff] = 0.0
    
    nonzero = dotp[~np.all(dotp == 0, axis=1)]

    return nonzero.T


def dihedral(p1,p2,p3,p4):
    """
    calculate dihedral angle given 4 position vectors
    Parameters
    ----------
    p1,p2,p3,p4: (3) numpy array 
        position vectors

    Returns
    -------
    dotp : float
        angle in radians, range (-pi,pi)
    
    """
    
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




