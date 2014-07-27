#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014 Sandro Bottaro (sbottaro@sissa.it)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as N
from scipy.spatial import distance


a = 5.0
b = 5.0
c = 3.0
scale = [1.0/a,1.0/b,1.0/c]
interactions = ['..','>>','<<','<>','><','WC','WW','WS','WH','HH','HS','HW','SS','SH','SW','XX']



def coord2lcs(atoms):
    
    # calculate center of mass
    origo = N.sum(atoms,axis=1)/3.0
    # CoM-C2 (x axis)
    x = atoms[:,0]-origo
    x_norm = 1./N.sqrt(N.sum(x*x,axis=1))
    x = x*N.transpose([x_norm,x_norm,x_norm])
    # CoM-C4/C6 
    c = atoms[:,1]-origo
    # z/y axis
    z = N.cross(x,c,axis=1)
    z_norm = 1./N.sqrt(N.sum(z*z,axis=1))
    z = z*N.transpose([z_norm,z_norm,z_norm])
    y = N.cross(z,x,axis=1)
    lcs = N.array([x.T,y.T,z.T]).T
    return lcs,origo


# standard version 

def lcs2mat_3d(lcs,origo,cutoff):

    ll = len(origo)
    # mat: x,y,z,d
    mat = N.zeros((ll,ll,3))

    # prune search first
    max_r  = N.sqrt(a*a*cutoff*cutoff)
    dmat = distance.pdist(origo)
    c_idx = (dmat<max_r).nonzero()[0]
    m_idx = N.array(N.triu_indices(ll,1)).T[c_idx] 

    # loop over pruned array
    for kk in xrange(len(c_idx)):
        ii = m_idx[kk][0]
        jj = m_idx[kk][1]
        local_cs1 = N.matrix(lcs[ii])
        local_cs2 = N.matrix(lcs[jj])
        p1 = origo[ii]
        p2 = origo[jj]
        diff = p2-p1

        R1 = diff*local_cs1
        R2 = -diff*local_cs2
        R1_scaled = N.array([scale[l]*R1[0,l] for l in range(3)])
        R2_scaled = N.array([scale[l]*R2[0,l] for l in range(3)])
        
        D1 = N.sqrt(sum(R1_scaled**2))
        D2 = N.sqrt(sum(R2_scaled**2))

        if(D1 < cutoff):
            for nn in range(3):
                mat[ii,jj,nn] = R1[0,nn]

        if(D2 < cutoff):
            for nn in range(3):
                mat[jj,ii,nn] = R2[0,nn]

    return mat

def lcs2mat_score(lcs,origo,cutoff):

    ll = len(origo)
    mat = []
    ids = []

    # prune search first
    max_r  = N.sqrt(a*a*cutoff*cutoff)
    dmat = distance.pdist(origo)
    c_idx = (dmat<max_r).nonzero()[0]
    m_idx = N.array(N.triu_indices(ll,1)).T[c_idx] 

    # loop over pruned array
    for kk in xrange(len(c_idx)):
        ii = m_idx[kk][0]
        jj = m_idx[kk][1]
        local_cs1 = N.matrix(lcs[ii])
        local_cs2 = N.matrix(lcs[jj])
        p1 = origo[ii]
        p2 = origo[jj]
        diff = p2-p1

        R1 = diff*local_cs1
        R2 = -diff*local_cs2
        R1_scaled = N.array([scale[l]*R1[0,l] for l in range(3)])
        R2_scaled = N.array([scale[l]*R2[0,l] for l in range(3)])
        
        D1 = N.sqrt(sum(R1_scaled**2))
        D2 = N.sqrt(sum(R2_scaled**2))

        if(D1 < cutoff):
            v = []
            for nn in range(3):
                v.append(R1[0,nn])
            mat.append(v)
            ids.append([ii,jj])

        if(D2 < cutoff):
            v = []
            for nn in range(3):
                v.append(R2[0,nn])
            mat.append(v)
            ids.append([jj,ii])

    return N.array(mat).T,ids


def lcs2mat_3d_scaled(lcs,origo,cutoff):

    ll = len(origo)
    # mat: x,y,z,d
    mat = N.zeros((ll,ll,3))

    # prune search first
    max_r  = N.sqrt(a*a*cutoff*cutoff)
    dmat = distance.pdist(origo)
    c_idx = (dmat<max_r).nonzero()[0]
    m_idx = N.array(N.triu_indices(ll,1)).T[c_idx] 

    # loop over pruned array
    for kk in xrange(len(c_idx)):
        ii = m_idx[kk][0]
        jj = m_idx[kk][1]
        local_cs1 = N.matrix(lcs[ii])
        local_cs2 = N.matrix(lcs[jj])
        p1 = origo[ii]
        p2 = origo[jj]
        diff = p2-p1

        R1 = diff*local_cs1
        R2 = -diff*local_cs2
        R1_scaled = N.array([scale[l]*R1[0,l] for l in range(3)])
        R2_scaled = N.array([scale[l]*R2[0,l] for l in range(3)])
        
        D1 = N.sqrt(sum(R1_scaled**2))
        D2 = N.sqrt(sum(R2_scaled**2))

        if(D1 < cutoff):
            for nn in range(3):
                mat[ii,jj,nn] = R1_scaled[nn]

        if(D2 < cutoff):
            for nn in range(3):
                mat[jj,ii,nn] = R2_scaled[nn]

    return mat

def lcs2mat_1d(lcs,origo,cutoff):

    ll = len(origo)
    # mat: x,y,z,d
    mat_sb = N.zeros((ll,ll))

    # prune search first
    max_r  = N.sqrt(a*a*cutoff*cutoff)
    dmat = distance.pdist(origo)
    c_idx = (dmat<max_r).nonzero()[0]
    m_idx = N.array(N.triu_indices(ll,1)).T[c_idx] 

    # loop over pruned array
    for kk in xrange(len(c_idx)):
        ii = m_idx[kk][0]
        jj = m_idx[kk][1]
        local_cs1 = N.matrix(lcs[ii])
        local_cs2 = N.matrix(lcs[jj])
        p1 = origo[ii]
        p2 = origo[jj]
        diff = p2-p1

        R1 = diff*local_cs1
        R2 = -diff*local_cs2
        R1_scaled = N.array([scale[l]*R1[0,l] for l in range(3)])
        R2_scaled = N.array([scale[l]*R2[0,l] for l in range(3)])
        
        D1 = N.sqrt(sum(R1_scaled**2))
        D2 = N.sqrt(sum(R2_scaled**2))

        if(D1 < cutoff):
            mat_sb[ii,jj] = cutoff - D1

        if(D2 < cutoff):
            mat_sb[jj,ii] = cutoff - D2

    return mat_sb



def lcs2mat_4d(lcs,origo,cutoff):

    ll = len(origo)
    mat_gb = N.zeros((ll,ll,4))

    # prune search first
    max_r  = N.sqrt(a*a*cutoff*cutoff)
    dmat = distance.pdist(origo)
    c_idx = (dmat<max_r).nonzero()[0]
    m_idx = N.array(N.triu_indices(ll,1)).T[c_idx] 

    # loop over pruned array
    for kk in xrange(len(c_idx)):
        ii = m_idx[kk][0]
        jj = m_idx[kk][1]
        local_cs1 = N.matrix(lcs[ii])
        local_cs2 = N.matrix(lcs[jj])
        p1 = origo[ii]
        p2 = origo[jj]
        diff = p2-p1

        R1 = diff*local_cs1
        R2 = -diff*local_cs2
        R1_scaled = N.array([scale[l]*R1[0,l] for l in range(3)])
        R2_scaled = N.array([scale[l]*R2[0,l] for l in range(3)])
        
        D1 = N.sqrt(sum(R1_scaled**2))
        D2 = N.sqrt(sum(R2_scaled**2))

        if(D1 < cutoff):
            D1_S = (D1*N.pi)/cutoff
            s = N.sin(D1_S)/D1_S
            for nn in range(3):
                mat_gb[ii,jj,nn] = s*R1_scaled[nn]
                    
            mat_gb[ii,jj,3] = 1.0+N.cos(D1_S)

        if(D2 < cutoff):
            D2_S = (D2*N.pi)/cutoff
            s = N.sin(D2_S)/D2_S
            for nn in range(3):
                mat_gb[jj,ii,nn] = s*R2_scaled[nn]
                    
            
            #    mat_gb[jj,ii,nn] = s*R2_scaled[nn]
            mat_gb[jj,ii,3] = 1.0+N.cos(D2_S)

    return mat_gb



def calc_dist_nd(m1,m2):

    
    assert(m1.shape==m2.shape)
    ll = m1.shape[0]
    diff = (m1-m2)**2
    norm = N.sum(diff,axis=2)
    #dd = (N.sum(norm,axis=1) + N.sum(norm,axis=0))/(2.0*ll)
    ermsd = N.sqrt(sum(sum(norm))/ll)
    return ermsd


def calc_dist_1d(m1,m2):

    assert(m1.shape==m2.shape)
    ll = m1.shape[0]
    diff = (m1-m2)**2
    #dd = (N.sum(diff,axis=1) + N.sum(diff,axis=0))/(2.0*ll)
    ermsd = N.sqrt(sum(sum(diff))/ll)
    return ermsd



# analyze mat file 
def analyze_mat(M,seq):
    

    # mean values and covariance matrix for wc-pair calculation
    # extracted from empirical distribution
    mean = [2.86,4.67,0.01]
    sigma = [[ 0.26, -0.14,0.0],\
                 [-0.14,0.13,0.0],\
                 [ 0.0, 0.0,  0.33 ]]
    # treshold values for base pair edges were obtained
    # from the angular distribution
    
    theta1 = 0.16
    theta2 = 2.0
    theta3 = -2.0


    def multivariate_pdf(r, mean, cov):
        dim  = r.shape[-1]
        dev  = r - mean
        maha = N.einsum('...k,...kl,...l->...', dev, N.linalg.pinv(cov), dev)
        return (2 * N.pi)**(-0.5 * dim) * N.linalg.det(cov)**(0.5) * N.exp(-0.5 * maha)

    def is_stack(v1,v2):
        
        rho1 = N.sqrt(v1[0]**2 + v1[1]**2)
        rho2 = N.sqrt(v2[0]**2 + v2[1]**2)

        if(rho1 > 5.0 and rho2 > 5.0):
            return 0
        if(v1[2] > 2.0 and v2[2] > 2.0):
            return interactions.index('><')
        if(v1[2] > 2.0 and v2[2] < -2.0):
            return interactions.index('>>')
        if(v1[2] < -2.0 and v2[2] > 2.0):
            return interactions.index('<<')
        if(v1[2] < -2.0 and v2[2] < -2.0):
            return interactions.index('<>')
        return 0

    def is_wc(v1,v2,r1,r2):
        if((r1 == 'A' and r2 == 'U') or  (r1 == 'U' and r2 == 'A') or \
           (r1 == 'C' and r2 == 'G') or  (r1 == 'G' and r2 == 'C')):
            rv1 = multivariate_pdf(v1,mean,sigma) 
            rv2 = multivariate_pdf(v2,mean,sigma) 
            if(rv1*rv2>1.0e-08):
                return interactions.index('WC')
        return 0

    def is_noncanonical(v1,v2):

        def edge(angle):
            if(angle > theta1 and angle <= theta2):
                return 'W'
            if(angle <= theta1 and angle > theta3):
                return 'S'
            return 'H'

        if(N.abs(v1[2]) > 2.0 or N.abs(v2[2]) > 2.0):
            return 0
        
        e1 = edge(N.arctan2(v1[1],v1[0]))
        e2 = edge(N.arctan2(v2[1],v2[0]))
        return interactions.index(e1+e2)

    ll = M.shape[0]
    int_mat = N.zeros((ll,ll))

    for i in xrange(ll):
        r1 = seq[i].split("_")[1]
        for j in xrange(i+1,ll):
            r2 = seq[j].split("_")[1]
            
            if(all(M[i,j]) == 0.0 or all(M[j,i]) == 0.0):
                continue

            # check if bases are stacked
            stack = is_stack(M[i,j],M[j,i])
            if(stack > 0):
                int_mat[i,j] = stack
                continue

            # check watson-crick pairing
            wc = is_wc(M[i,j],M[j,i],r1,r2) 
            if(wc >0):
                int_mat[i,j] = wc
                continue

            # check other pairing            
            nwc = is_noncanonical(M[i,j],M[j,i])    
            if(nwc >0):
                int_mat[i,j] = nwc
                continue

            # give generic
            int_mat[i,j] = interactions.index("XX")


    return int_mat


