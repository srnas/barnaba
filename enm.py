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


import reader as reader
import numpy as np
from scipy.linalg import eigh
from scipy.spatial import distance





def enm(args):
    TOL=0.000001 #Tollerance to identify zero-eigenvalues modes, 
               #corresponding to translational and rotational degrees of freedom.
               #Set it to zero or a negative value to print all the eigenvalues/eigenvectors


#    print "# Calculating ENM angles..."

    files = args.files

#    fh = open(args.name,'w')
#    fh.write("# This is a baRNAba run.\n")
#    for k in args.__dict__:
#        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
#        fh.write(s)

    if(args.type=="S"):
        atoms_req = ["C1'"]
    if(args.type=="B"):
        atoms_req = ["C2"]
    if(args.type=="P"):
        atoms_req = ["P"]
    if(args.type=="SBP"):
        atoms_req = ["P","C1'","C2"]
    if(args.type=="AA"):
        atoms_req = reader.Names.rna_atoms

    # loop over files
    for i in xrange(0,len(files)):
        
        cur_pdb = reader.Pdb(files[i],base_only=False)

        # loop over models
        for j in xrange(len(cur_pdb.models)):
            
            beads = []
            #find beads
            for k in xrange(len(cur_pdb.models[j].sequence)):
                resi = cur_pdb.models[j].residues[k]
                for atom_type in atoms_req:
                    cc = resi.get_atomobject(atom_type)
                    if(cc!=None):
                        beads.append(cc)
            print "# Read ", len(beads), "coordinates"

            coords=[]
            for b in beads:
                coords.append(b.coords)
            coords = np.array(coords)
            # build distance matrix
            dmat = distance.pdist(coords)
            ll = len(coords)
            c_idx = (dmat<args.cutoff).nonzero()[0]
            m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 

            # Build Interaction Matrix
            mat = np.zeros((3*ll,3*ll))
            k_elast=1.
            roba_i=[0,0,0,1,1,2]
            roba_j=[0,1,2,1,2,2]
            for kk in xrange(len(c_idx)):
                ii = m_idx[kk][0]
                jj = m_idx[kk][1]
                d=[]
                for mu in range(0,3):
                    d.append(coords[ii,mu]-coords[jj,mu])
                cacca=k_elast/dmat[c_idx[kk]]**2
                iii=3*ii
                jjj=3*jj
                for mu in [0,1,2]:                                        
                    temp=(d[mu]**2)*cacca
                    mat[iii+mu,iii+mu]+=temp # diagonale
                    mat[jjj+mu,iii+mu]+=-temp
                    mat[jjj+mu,jjj+mu]+=temp  # diagonale                
                for k in [1,2,4]:                                        
                    mu=roba_i[k]
                    nu=roba_j[k]
                    temp=d[mu]*d[nu]*cacca
                    mat[jjj+mu,iii+nu]+=-temp
                    mat[iii+nu,iii+mu]+=temp # diagonal
                    mat[jjj+nu,iii+mu]+=-temp
                    mat[jjj+nu,jjj+mu]+=temp  # diagonal               

            print '# Diagonalization'
            eigens=eigh(mat,lower=True)

            # fix phase to make tests reproducible
            for i in xrange(len(eigens[0])):
                m=0.0
                for k in xrange(len(eigens[1][:,i])):
                    if(eigens[1][k,i]>TOL):
                        m=1.0
                        break
                    if(eigens[1][k,i]<-TOL):
                        m=-1.0
                        break
                eigens[1][:,i]*=m

            print '# Writing to files'
            feval=open(args.name+".eval.dat",'w')
            for i in xrange(len(eigens[0])):
                if (not args.zmodes) and eigens[0][i]<TOL:
                    stri = "%5d %.6e \n" % (i,0.0)
                else:
                    stri = "%5d %.6e \n" % (i,eigens[0][i])
                feval.write(stri)
            feval.close()
            
            MAXVEC=args.ntop
            nvec=0
            for k in xrange(len(eigens[0])):
                if nvec>MAXVEC:
                    break
                if (not args.zmodes) and eigens[0][k]<TOL:
                    continue
                fevec=open(args.name+".evec"+str(k).zfill(len(str(args.ntop+6)))+".dat",'w')
                fevec.write("# Eigenvector number "+str(k)+"\n")
                fevec.write("# beads index & x-component & y-component & z-component \n")
                stri = ""
                for i in xrange(len(eigens[1][:,k])/3):                    
                    stri += "%5d %10.6f %10.6f %10.6f \n" % (i,eigens[1][3*i+0,k],eigens[1][3*i+1,k],eigens[1][3*i+2,k])
                fevec.write(stri) 
                fevec.close()
                nvec+=1

            # Print distance fluctuations between consecutive C2 atoms
            print "# Computing distance fluctuations between consecutive C2 atoms"
            i=0
            C2_indeces=[]
            for b in beads:
                if b.atom_t=="C2":
                    C2_indeces.append(i)
                i+=1
            print "# There are",len(C2_indeces),"C2 atoms"
            fC2=open(args.name+".distC2.dat",'w')
            fC2.write("# fluctuations of the distances between consecutive C2 atoms in the ENM model \n")
            fC2.write("# bead index i & bead index j & sigma^2 i--i+1 \n")
            for n in xrange(len(C2_indeces)-1):
                i=C2_indeces[n]
                j=C2_indeces[n+1]

                d=[]
                for mu in range(0,3):
                    d.append(coords[i,mu]-coords[j,mu])
                d=np.array(d)
                d/=np.sqrt(sum(d*d))
                sigma=0.0
                for  mu in range(0,3):
                    for  nu in range(0,3):
                        C_ii_munu=0
                        C_jj_munu=0
                        C_ij_munu=0
                        C_ji_munu=0
                        for k in xrange(len(eigens[0])):
                            if eigens[0][k]>TOL:
                                lambd=1./eigens[0][k]
                                v_i_mu=eigens[1][3*i+mu,k]
                                v_i_nu=eigens[1][3*i+nu,k]
                                v_j_mu=eigens[1][3*j+mu,k]
                                v_j_nu=eigens[1][3*j+nu,k]
                                C_ii_munu+=lambd*v_i_mu*v_i_nu
                                C_jj_munu+=lambd*v_j_mu*v_j_nu
                                C_ij_munu+=lambd*v_i_mu*v_j_nu
                                C_ji_munu+=lambd*v_j_mu*v_i_nu
                        sigma+=d[mu]*d[nu]*(C_ii_munu+C_jj_munu-C_ij_munu-C_ji_munu)
                    
                stri = "%5d %5d %10.6f \n" % (i,j,sigma)
                fC2.write(stri)
            fC2.close()

#    fh.close()
    return 0
