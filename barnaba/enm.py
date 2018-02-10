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
from scipy.spatial import distance
import mdtraj as md
from . import definitions


class Enm:

    def __init__(self,pdb,sele_atoms,cutoff=8.0):

        cur_pdb = md.load_pdb(pdb)
        topology = cur_pdb.topology
        if(sele_atoms=="AA"):
            idxs = [atom.index for atom in topology.atoms if ("H" not in atom.name)]
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

        # construct matrix
        mat = np.zeros((ll,ll,3,3))
        for kk in range(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            mat[ii,jj] = -k_elast[kk]*np.outer(diff[kk],diff[kk])

        # now fill the 3n x 3n matrix (there might be a more pythonic way. but this works)
        mat1 = np.zeros((3*ll,3*ll))
        for k1 in range(ll):
            # diagonal elements
            diag = np.sum(mat[k1,:],axis=0) + np.sum(mat[:,k1],axis=0)
            for i1 in range(3):
                for i2 in range(i1,3):
                    mat1[3*k1+i1,3*k1+i2] =  -diag[i1,i2]
            # off diagonal
            for k2 in range(k1+1,ll):
                for i1 in range(3):
                    for i2 in range(3):
                        mat1[3*k1+i1,3*k2+i2] =  mat[k1,k2,i1,i2]
                
        # diagonalise
        e_val,e_vec=eigh(mat1.T,lower=True)
        self.e_val = e_val
        self.e_vec = e_vec
        self.coords = coords
        #self.idx_c2 = [cur_pdb.topology.atom(idxs[x]).index for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")]
        self.idx_c2 = [x for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")]
        self.seq_c2 = [str(cur_pdb.topology.atom(idxs[x]).residue) for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")]


    def get_eval(self):
        return self.e_val

    def get_evec(self):
        return self.e_vec

    def c2_fluctuations(self):


        if(len(self.idx_c2)==0):
            print("# no C2 atoms in PDB")
            exit(1)

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
            # check eigenvalue to be nonzero
            assert(self.e_val[i] > definitions.tol)
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
