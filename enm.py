#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014 Sandro Bottaro (sbottaro@sissa.it) and Giovanni Pinamonti (giopina@sissa.it)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
from scipy.linalg import eigh
from scipy.spatial import distance
import mdtraj as md
import definitions



def enm(args):


    # define atoms
    if(args.type=="S"):
        atoms_req = "(name \"C1'\" and (resname " + " or resname ".join(definitions.rna) + "))"
    if(args.type=="B"):
        atoms_req = "(name C2 and (resname " + " or resname ".join(definitions.rna) + "))"
    if(args.type=="P"):
        atoms_req = "(name P and (resname " + " or resname ".join(definitions.rna) + "))"
    if(args.type=="SBP"):
        atoms_req = "(( name P or name \"C1'\" or name C2 ) and (resname " + " or resname ".join(definitions.rna) + " ))"
    if(args.type=="AA"):
        atoms_req = "( (name " + " or name ".join(definitions.heavy_atoms) + ") and (resname " + " or resname ".join(definitions.rna) + " ))"
        
    if(args.protein):
        atoms_req += "or name CA"
    cur_pdb = md.load_pdb(args.pdbs)
    idxs = cur_pdb.topology.select(atoms_req)

    coords = cur_pdb.xyz[0,idxs]


    print "# Read ", coords.shape, "coordinates"
    # build distance matrix
    dmat = distance.pdist(coords)
    ll = len(coords)

    # find where distance is shorter than cutoff
    c_idx = (dmat<args.cutoff).nonzero()[0]
    m_idx = np.array(np.triu_indices(ll,1)).T[c_idx]
    # k = gamma/d^2
    k_elast=1./dmat[c_idx]**2

    # difference
    diff = [coords[ii]-coords[jj] for ii,jj in m_idx]

    # construct matrix
    mat = np.zeros((ll,ll,3,3))
    for kk in xrange(len(c_idx)):
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

    # write to file - eigenvalues
    fh = open(args.name,'a')
    fh.write("# Eigenvalues \n")
    for i in xrange(len(e_val)):
        vv = e_val[i]
        if(vv<definitions.tol): vv = 0.0
        stri = "%5d %.6e \n" % (i,vv)
        fh.write(stri)
    fh.close()
    
    # write eigenvectors (skip first 6)
    for i in xrange(6,args.ntop+6):
        # check eigenvalue to be nonzero
        assert(e_val[i] > definitions.tol)
        fh= open(args.name + "." + str(i).zfill(2),'w')
        fh.write("# eigenvector " + str(i) + "\n")
        fh.write("# beads index & x-component & y-component & z-component \n")
        stri = ""
        ee = e_vec[:,i]
        # check phase (this makes tests reproducible)
        if(ee[0]<definitions.tol): ee*= -1.0
        for k in xrange(ll):                    
            stri += "%5d %10.6f %10.6f %10.6f \n" % (k,ee[3*k],ee[3*k+1],ee[3*k+2])
        fh.write(stri) 
        fh.close()

    
    if(args.type == "P" or args.type == "S"):
        return 0
    
    # C2-C2 fluctuations
    C2_indeces = [x for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")]

    # assume that C2 are all consecutive - maybe a check on the chain would be useful?
    if(len(C2_indeces)==0):
        print "# no C2 atoms in PDB"
        return 0

    # get coordinates
    fh = open(args.name + ".SHAPE","w")
    
    for n in range(len(C2_indeces)-1):
        i = 3*C2_indeces[n]
        j = 3*(C2_indeces[n+1])
        diff = coords[C2_indeces[n]]-coords[C2_indeces[n+1]] 
        diff /= np.sqrt(np.sum(diff**2))
        
        v_i = [e_vec[i:i+3,k] for k in xrange(6,len(e_val))]
        v_j = [e_vec[j:j+3,k] for k in xrange(6,len(e_val))]

        top = len(e_val)-6
        
        c_ii = np.array([np.outer(v_i[k],v_i[k])/e_val[k+6] for k in xrange(top)])
        c_jj = np.array([np.outer(v_j[k],v_j[k])/e_val[k+6] for k in xrange(top)])
        c_ij = np.array([np.outer(v_i[k],v_j[k])/e_val[k+6] for k in xrange(top)])
        c_ji = np.array([np.outer(v_j[k],v_i[k])/e_val[k+6] for k in xrange(top)])

        # sum contributions from all eigenvectors
        tensor = np.sum([(c_ii[k]+c_jj[k] - c_ij[k] - c_ji[k]) for k in xrange(top)],axis=0)
        sigma = np.dot(diff,np.dot(tensor,diff))
                                        
        fh.write("%5i %5d %10.6f \n" % (i,j,sigma))
        
    fh.close()
    
    return 0
