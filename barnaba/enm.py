#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2017 Sandro Bottaro (sandro.bottaro@bio.ku.dk) and Giovanni Pinamonti (GP)

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

# Make sure that range returns an iterator also in python2 (using future module)
from builtins import range

import numpy as np
from scipy.linalg import eigh
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import eigsh
from scipy.spatial import distance
import mdtraj as md
import sys
from . import definitions


class Enm:
    '''Creates an elastic network model and diagonalizes the interaction matrix to find the 
    principal modes.
    See Pinamonti et al., NAR 2015 for a more detailed description.
    ------------
     parameters
    ------------
    pdb        : mdtraj trajectory object (TODO: what happens with multiple frames?)
    sele_atoms : atoms to use as beads (default=["C1\'","C2","P","CA","CB"])
    cutoff     : cutoff radius in nm (default=0.9)
    sparse     : whether or not to use sparse matrices in the diagonalization (default=False)
    ntop       : number of eigenvectors to print, excluding the ones corresponding to null eigenvalues (default=10)
    '''
    def __init__(self,pdb,sele_atoms=["C1\'","C2","P","CA","CB"],cutoff=0.9,sparse=False,ntop=10):
        self.sparse=sparse
        cur_pdb = md.load_pdb(pdb)
        topology = cur_pdb.topology
        if(sele_atoms=="AA"):
            idxs = [atom.index for atom in topology.atoms if (atom.name[0] != "H" and not ( atom.name[0].isdigit() and atom.name[1] == "H") )]
        else:
            idxs = [atom.index for atom in topology.atoms if (atom.name in sele_atoms)]
        if(len(idxs)==0):
            print("# Error. no atoms found.")
            exit(1)
            
        # define atoms
        coords = cur_pdb.xyz[0,idxs]

        print("# Read ", coords.shape, "coordinates")
        # build distance matrix
        dmat = distance.pdist(coords)

        ll = len(coords)

        # find where distance is shorter than cutoff
        c_idx = (dmat<cutoff).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx]
        # k = gamma/d^2
        k_elast=1./dmat[c_idx]**2
        
        # difference
        diff = [coords[ii]-coords[jj] for ii,jj in m_idx]
        if self.sparse:
            # construct matrix
            ele_up=np.zeros(len(c_idx)*9)
            idx_up=np.zeros((2,len(c_idx)*9))
            ele_diag=np.zeros(ll*9)
            idx_diag=np.zeros((2,ll*9))
            kkk=0
            for i in range(ll):
                for mu in range(3):
                    for nu in range(3):
                        idx_diag[:,kkk]=(3*i+mu,3*i+nu)
                        kkk+=1
            kkk=0
            for k in range(len(c_idx)):
                i = m_idx[k][0]
                j = m_idx[k][1]
                for mu in range(3):
                    for nu in range(3):
                        temp=k_elast[k]*diff[k][mu]*diff[k][nu]
                        # filling off-diagonal elements
                        ele_up[kkk]=-temp
                        idx_up[:,kkk]=(3*i+mu,3*j+nu)
                        kkk+=1
                        # filling diagonal elements
                        ele_diag[9*i+3*mu+nu]+=temp
                        ele_diag[9*j+3*mu+nu]+=temp
            idx_down=np.array((idx_up[1],idx_up[0]))
            idx_tot=np.concatenate([idx_up,idx_down,idx_diag],axis=1)
            ele_tot=np.concatenate([ele_up,ele_up,ele_diag])
            mat=sp.csc_matrix((ele_tot,idx_tot))
            self.mat=mat # store the interaction matrix for future C2-C2 calculations
            # diagonalise
            print('# Using sparse matrix diagonalization')
            e_val,e_vec=eigsh(mat, k=ntop+6,sigma=definitions.tol)
        else:
            # construct matrix
            mat = np.zeros((3*ll,3*ll))
            for kk in range(len(c_idx)):
                ii = m_idx[kk][0]
                jj = m_idx[kk][1]
                mat[3*ii:3*ii+3,3*jj:3*jj+3] = -k_elast[kk]*np.outer(diff[kk],diff[kk])
                mat[3*ii:3*ii+3,3*ii:3*ii+3] += k_elast[kk]*np.outer(diff[kk],diff[kk])
                mat[3*jj:3*jj+3,3*jj:3*jj+3] += k_elast[kk]*np.outer(diff[kk],diff[kk])
                
            # diagonalise
            e_val,e_vec=eigh(mat.T,lower=True)
        ### check here:
        ### 2) MAXVEC has to be obtained from args.ntop
        ###    Done. Do I want to print the 0 modes or not?
        ### 4) EIGSH IS GIVIN EXTRA ZERO-MODES!
        ###    Sigma has to be >zero to avoid extra null modes to pop out
        ###    if sigma > 10x smallest eval => wrong results
        ###    I set sigma=tol. This should work if tol makes sense
        self.check_null_modes(e_val,ntop)

        
        self.e_val = e_val ### GP Is there a particular reason for not doing this before
        self.e_vec = e_vec
        self.coords = coords
        #self.idx_c2 = [cur_pdb.topology.atom(idxs[x]).index for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")]
        # get C2 indexes for future C2-C2 fluctuations
        self.idx_c2 = np.array([x for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")])
        self.seq_c2 = [str(cur_pdb.topology.atom(idxs[x]).residue) for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")]


    def check_null_modes(self,e_val,ntop):
        N_NULL=6
        for i in range(6,ntop+6):
            if(e_val[i] < definitions.tol):
                N_NULL+=1
        if N_NULL>6:
            sys.stderr.write("WARNING: there are %d null modes. \
            Normally there should be only 6 corresponding to rotational and translational invariance.\
            This can lead to unpredictable results." % N_NULL)
        
    def get_eval(self):
        return self.e_val

    def get_evec(self):
        return self.e_vec

    def c2_fluctuations(self):

        # check if there are C2 atoms in the structure (needed to compute C2-C2 fluctuations...)
        if(len(self.idx_c2)==0):
            print("# no C2 atoms in PDB: can't compute C2-C2 fluctuations")
            exit(1)

        if self.sparse:
            # To deal with big matrices the effective interaction matrix
            # of the C2 beads is computed (see Pinamonti et al. NAR 2015)
            ll=len(self.coords)
            if len(self.idx_c2)*3==ll :
                #skip effective interaction computation if there are only C2 beads
                M_new=mat
            else:
                # submatrices
                idx_a=np.sort(np.concatenate([3*self.idx_c2+i for i in range(3)]))
                M_a=self.mat[np.ix_(idx_a,idx_a)]
                idx_b=np.delete(np.arange(3*ll),idx_a)
                M_b=self.mat[np.ix_(idx_b,idx_b)]
                W=self.mat[np.ix_(idx_a,idx_b)]
                
                X=spsolve(M_b.tocsc(),W.T.tocsc())
                ### X is not a sparse matrix... Damn...
                # effective interaction matrix between C2 beads
                M_new=(M_a-sp.csc_matrix.dot(W.tocsc(),X)).toarray()
                eval_new,evec_new=eigh(M_new)
            # Computing the fluctuations
            sigma=[]
            for n in range(len(self.idx_c2)-1):
                i = 3*n
                j = 3*(n+1)
                diff = self.coords[self.idx_c2[n]]-self.coords[self.idx_c2[n+1]] 
                diff /= np.sqrt(np.sum(diff**2))
                
                v_i = [evec_new[i:i+3,k] for k in range(6,len(eval_new))]
                v_j = [evec_new[j:j+3,k] for k in range(6,len(eval_new))]
                
                top = len(eval_new)-6
                
                c_ii = np.array([np.outer(v_i[k],v_i[k])/eval_new[k+6] for k in range(top)])
                c_jj = np.array([np.outer(v_j[k],v_j[k])/eval_new[k+6] for k in range(top)])
                c_ij = np.array([np.outer(v_i[k],v_j[k])/eval_new[k+6] for k in range(top)])
                c_ji = np.array([np.outer(v_j[k],v_i[k])/eval_new[k+6] for k in range(top)])
                
                # sum contributions from all eigenvectors
                tensor = np.sum([(c_ii[k]+c_jj[k] - c_ij[k] - c_ji[k]) for k in range(top)],axis=0)
                sigma.append(np.dot(diff,np.dot(tensor,diff)))

        else:

            sigma = []
            for n in range(len(self.idx_c2)-1):
                i = 3*self.idx_c2[n]
                j = 3*(self.idx_c2[n+1])
                diff = self.coords[self.idx_c2[n]]-self.coords[self.idx_c2[n+1]] 
                diff /= np.sqrt(np.sum(diff**2))
                
                v_i = [self.e_vec[i:i+3,k] for k in range(6,len(self.e_val))]
                v_j = [self.e_vec[j:j+3,k] for k in range(6,len(self.e_val))]
                
                top = len(self.e_val)-6
                
                c_ii = np.array([np.outer(v_i[k],v_i[k])/self.e_val[k+6] for k in range(top)])
                c_jj = np.array([np.outer(v_j[k],v_j[k])/self.e_val[k+6] for k in range(top)])
                c_ij = np.array([np.outer(v_i[k],v_j[k])/self.e_val[k+6] for k in range(top)])
                c_ji = np.array([np.outer(v_j[k],v_i[k])/self.e_val[k+6] for k in range(top)])
            
                # sum contributions from all eigenvectors
                tensor = np.sum([(c_ii[k]+c_jj[k] - c_ij[k] - c_ji[k]) for k in range(top)],axis=0)
                sigma.append(np.dot(diff,np.dot(tensor,diff)))
        return sigma, self.seq_c2


    def print_eval(self):
        stri = "# Eigenvalues \n"
        for i in range(len(self.e_val)):
            vv = self.e_val[i]
            if(vv<definitions.tol): vv = 0.0
            stri += "%5d %.6e \n" % (i,vv)
        return stri
    
    def print_evec(self,ntop):

        stri = ""
        # write eigenvectors (skip first 6)
        for i in range(6,ntop+6):
            if(self.e_val[i] < definitions.tol):
               stri += "# skipping eigenvector %d (eigenvalue = %f )\n" % (i,self.e_val[i])
               continue
            stri += "# eigenvector %d \n" % (i) 
            stri += "# beads index & x-component & y-component & z-component \n"
            
            ee = self.e_vec[:,i]
            # check phase (this makes tests reproducible)
            if(ee[0]<definitions.tol): ee*= -1.0
            for k in range(len(self.coords)):                    
                stri += "%5d %10.6f %10.6f %10.6f \n" % (k,ee[3*k],ee[3*k+1],ee[3*k+2])
            stri += "\n"
            stri += "\n"
        return stri
