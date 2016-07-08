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
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import eigsh
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
    print 'dim',ll,3*ll ###
    # find where distance is shorter than cutoff
    c_idx = (dmat<args.cutoff).nonzero()[0]
    m_idx = np.array(np.triu_indices(ll,1)).T[c_idx]
    # k = gamma/d^2
    k_elast=1./dmat[c_idx]**2

    # difference
    diff = [coords[ii]-coords[jj] for ii,jj in m_idx]

    MAXVEC=args.ntop+6
    import time
    t=time.time()

    if args.sparse:
        # construct matrix
        mat = sp.lil_matrix((ll*3,ll*3))
        #mat = dok_matrix((ll*3,ll*3)) ### slower
        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            ###
            for mu in range(3):
                for nu in range(3):
                    temp=k_elast[kk]*diff[kk][mu]*diff[kk][nu]
                    mat[3*ii+mu,3*jj+nu]+=-temp
                    mat[3*jj+mu,3*ii+nu]+=-temp
                    mat[3*ii+mu,3*ii+nu]+=temp
                    mat[3*jj+mu,3*jj+nu]+=temp

        print '# Using sparse matrix diagonalization'
        e_val,e_vec=eigsh(mat, k=MAXVEC,sigma=definitions.tol)
    else:
        mat = np.zeros((3*ll,3*ll))
        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            for mu in range(3):
                for nu in range(3):
                    temp=k_elast[kk]*diff[kk][mu]*diff[kk][nu]
                    mat[3*ii+mu,3*jj+nu]+=-temp
                    #mat[3*jj+mu,3*ii+nu]+=-temp ### no need to fill the lower triangle
                    mat[3*ii+mu,3*ii+nu]+=temp
                    mat[3*jj+mu,3*jj+nu]+=temp

        e_val,e_vec=eigh(mat.T,lower=True)
    print 'time',time.time()-t ###
### check here:
### 2) MAXVEC has to be obtained from args.ntop
###    Done. Do I want to print the 0 modes or not?
### 4) EIGSH IS GIVIN EXTRA ZERO-MODES!
###    Sigma has to be >zero to avoid extra null modes to pop out
###    if sigma > 10x smallest eval => wrong results
###    I set sigma=tol. This should work if tol makes sense

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
        ### we should put some message that explain what happens and tells to increase the cutoff here
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

    if not args.comp_shape:
        return 0
    if(args.type == "P" or args.type == "S"):
        return 0

    print '# computing C2-C2 fluctuations'
    # C2-C2 fluctuations
    C2_indeces = np.array([x for x in range(len(idxs)) if(cur_pdb.topology.atom(idxs[x]).name=="C2")])

    # assume that C2 are all consecutive - maybe a check on the chain would be useful?
    if(len(C2_indeces)==0):
        print "# no C2 atoms in PDB"
        return 0

    tshape=time.time()
    # get coordinates
    fh = open(args.name + ".SHAPE","w")
    if args.sparse:
        if args.type== "C2" :
            #skip effective interaction computation if there are only C2 beads
            ### maybe should test whether there are beads different from C2 or not
            M_new=mat
        else:
            # submatrices
            idx_a=np.sort(np.concatenate([3*C2_indeces+i for i in range(3)]))        
            M_a=mat[np.ix_(idx_a,idx_a)]
            idx_b=np.delete(np.arange(3*ll),idx_a)
            M_b=mat[np.ix_(idx_b,idx_b)]
            W=mat[np.ix_(idx_a,idx_b)]
            t=time.time()
            
            X=spsolve(M_b.tocsc(),W.T.tocsc())
            print '### time to spsolve =',time.time()-t
            ### X is not a sparse matrix... Damn...
            # effective interaction matrix between C2 beads
            M_new=(M_a-sp.csc_matrix.dot(W.tocsc(),X)).toarray()
            t=time.time()
            eval_new,evec_new=eigh(M_new)
            print '### time to eigh(M_new) =',time.time()-t

        for n in range(len(C2_indeces)-1):
            i = 3*n
            j = 3*(n+1)
            diff = coords[C2_indeces[n]]-coords[C2_indeces[n+1]] 
            diff /= np.sqrt(np.sum(diff**2))
            
            v_i = [evec_new[i:i+3,k] for k in xrange(6,len(eval_new))]
            v_j = [evec_new[j:j+3,k] for k in xrange(6,len(eval_new))]

            top = len(eval_new)-6
        
            c_ii = np.array([np.outer(v_i[k],v_i[k])/eval_new[k+6] for k in xrange(top)])
            c_jj = np.array([np.outer(v_j[k],v_j[k])/eval_new[k+6] for k in xrange(top)])
            c_ij = np.array([np.outer(v_i[k],v_j[k])/eval_new[k+6] for k in xrange(top)])
            c_ji = np.array([np.outer(v_j[k],v_i[k])/eval_new[k+6] for k in xrange(top)])

            # sum contributions from all eigenvectors
            tensor = np.sum([(c_ii[k]+c_jj[k] - c_ij[k] - c_ji[k]) for k in xrange(top)],axis=0)
            sigma = np.dot(diff,np.dot(tensor,diff))
                                        
            fh.write("%5i %5d %10.6f \n" % (i,j,sigma))
    
    else:
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
    print 'shape',time.time()-tshape

    return 0
