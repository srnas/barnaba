
import numpy as np
from scipy.spatial import distance
import residue as rd
import tools
import definitions
import sys


class Model:

    def __init__(self,data,set_lcs=True):

        ll = len(data)
        if(ll==0):
            sys.stderr.write("# Fatal error. No valid residues in pdb file \n")
            sys.exit(1)

        # array with all coordinates
        self.coords = np.array([data[i][5:8] for i in range(ll)])

        # list with sequence
        self.sequence = []
        self.sequence_id = []
        self.chain = []
        self.residues = []
        resi_tmp = []

        for i in xrange(ll):
            if(i==ll-1):
                next_resid = "XXX"
            else:
                next_resid = data[i+1][9]
            resi_tmp.append(data[i])
            if(next_resid != data[i][9]):
                self.residues.append(rd.Residue(resi_tmp,i))
                self.sequence.append(data[i][8])
                self.sequence_id.append(data[i][9])
                resi_tmp = []

        # set indeces for lcs construction
        if(set_lcs):
            self.set_lcs_idx()

        
    def __getitem__(self, ind):
        return self.residues[ind]

    def set_coords(self,data):
        self.coords = data

    # get ordered indeces for fast array accession
    def set_lcs_idx(self):
        
        CA_idx = []
        CB_idx = []
        CC_idx = []
        for i in xrange(len(self.sequence)):
            i1 = self.residues[i].get_idx("C2") 
            assert i1 != None, "# missing C2 in %s \n" % (self.sequence_id[i])
            i2 = self.residues[i].get_idx("C4") 
            assert i2 != None, "# missing C4 in %s \n" % (self.sequence_id[i])
            i3 = self.residues[i].get_idx("C6") 
            assert i3 != None, "# missing C6 in %s \n" % (self.sequence_id[i])

            CA_idx.append(i1)
            mytype = self.residues[i].res_mytype
            if(mytype in definitions.pyr):
                CB_idx.append(i2)
                CC_idx.append(i3)
            if(mytype in definitions.pur):
                CB_idx.append(i3)
                CC_idx.append(i2)

        self.lcs_idx =  np.array([CA_idx,CB_idx, CC_idx])


    # Find index of backbone atoms 
    def set_bb_index(self):

        # handle exceptions
        self.missing = [False]*6*len(self.sequence)
        # first and last are always missing
        self.missing[0] = True
        self.missing[-1] = True
        self.missing[-2] = True
        self.chi_missing = [False]*len(self.sequence)
        
        bb_idx_tmp = []
        chi_idx = []
        for i in xrange(len(self.sequence)):
            residue = self.residues[i]
            atoms = [residue.get_idx(at_type) for at_type in definitions.rna_backbone]
            for jj in xrange(len(atoms)):
                # flag when meissing
                if(atoms[jj]==None):
                    nos = [x for x in range(6*i+jj-3,6*i+jj+1) if x>-1]
                    for el in nos:
                        self.missing[el+1] = True
            bb_idx_tmp.extend(atoms)
            
            # now for chi angle
            atoms = [residue.get_idx("O4'"),residue.get_idx("C1'")]
            if(self.sequence[i] in definitions.pur):
                atoms.append(residue.get_idx("N9"))
                atoms.append(residue.get_idx("C4"))
            if(self.sequence[i] in definitions.pyr):
                atoms.append(residue.get_idx("N1"))
                atoms.append(residue.get_idx("C2"))
            if(None not in atoms):
                chi_idx.append(atoms)
            else:
                self.missing[i] = True
                

        # handle missing atoms
        bb_idx = []
        for ii in xrange(len(bb_idx_tmp)-3):
            idxs = bb_idx_tmp[ii:ii+4]
            if(None not in idxs):
                bb_idx.append(idxs)
        bb_angles = self.bb_idx = np.array(bb_idx)
        chi_angles = self.chi_idx = np.array(chi_idx)

        
    # calculate chi and backbone torsion angles
    def calc_bb_torsion(self):
        atom_pos_bb = self.coords[self.bb_idx]
        atom_pos_chi = self.coords[self.chi_idx]
        bb_angles = np.degrees(tools.dihedral(atom_pos_bb))
        chi_angles = np.degrees(tools.dihedral(atom_pos_chi))
        return bb_angles,chi_angles


    def set_j3_index(self):

        j3_idx = []
        lj = len(definitions.j3)
        self.j3_missing = [False]*lj*len(self.sequence)
        for i in xrange(len(self.sequence)):
            residue = self.residues[i]
            for j in xrange(len(definitions.j3)):
                
                atoms = [residue.get_idx(at_type) for at_type in definitions.j3[j][1]]
                if(definitions.j3[j][0]=="H3P"):
                    if(i!=len(self.sequence)-1):
                        atoms[-1] = self.residues[i+1].get_idx(definitions.j3[j][1][-1])
                    else:
                        # this is at position ll
                        atoms[-1] = None
                if(None in atoms):
                    self.j3_missing[i*lj+j] = True
                else:
                    j3_idx.append(atoms)
        self.j3_idx = np.array(j3_idx)

    def calc_j3(self):
        atom_pos = self.coords[self.j3_idx]
        angles = tools.dihedral(atom_pos)
        return np.degrees(angles)
        
    def set_pucker_index(self):
        
        pucker_idx = []
        
        self.pucker_missing = [False]*len(self.sequence)
        for i in xrange(len(self.sequence)):
            residue = self.residues[i]
            atoms = [residue.get_idx(at_type) for at_type in definitions.rna_pucker]

            # handle missing
            index = [ [atoms[k%5] for k in range(j,j+4)] for j in range(5)]
            findex = [item for sublist in index for item in sublist]
            if(None in findex):
                self.pucker_missing[i] = True
            else:
                pucker_idx.extend(index)
        self.pucker_idx = np.array(pucker_idx)

    def calc_pucker(self):
        atom_pos = self.coords[self.pucker_idx]
        angles = tools.dihedral(atom_pos)
        angles = angles.reshape(-1,5)
        
        # calculate pseudorotation parameters
        x1 = angles[:,4]+angles[:,1]-angles[:,3]-angles[:,0]

        #Pconst = 2.0*(np.sin(np.pi/5.) + np.sin(np.pi/2.5))
        x2 = 3.0776835*angles[:,2]
        
        # phase
        p0 = np.arctan2(x1,x2)
        
        p0[np.where( p0 < 0.0 )] += 2.0*np.pi
        
        # amplitude
        tm = angles[:,2]/np.cos(p0)

        #print angles.shape,p0.shape,tm.shape
        angles = np.concatenate((angles,tm[:,np.newaxis]),axis=1)
        angles = np.concatenate((angles,p0[:,np.newaxis]),axis=1)

        return np.degrees(angles)

    def set_h_idx(self):
        self.h_labels = []
        h_index = []
        for i in xrange(len(self.sequence)):
            residue = self.residues[i]
            for atom in residue.atom_types:
                if("H" in atom):
                    self.h_labels.append(self.sequence_id[i] + "_" +atom)
                    h_index.append(residue.get_idx(atom))
        assert len(self.h_labels)>0, "# Fatal error: your file does not contain Hydrogen atoms!"
        self.h_index = np.array(h_index)
        
    def calc_pairwise_h(self):
        atom_pos = self.coords[self.h_index]
        dmat = np.power(distance.pdist(atom_pos),-6)
        return dmat

        
    def get_lcs(self,ii):

        coords1 = self.coords[self.lcs_idx[0][ii]]
        coords2 = self.coords[self.lcs_idx[1][ii]]
        coords3 = self.coords[self.lcs_idx[2][ii]]
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

    def get_3dmat(self,cutoff,ii):

        lcs, origo = self.get_lcs(ii)

        # prune search first
        max_r  = np.max(definitions.f_factors)*cutoff
        dmat = distance.squareform(distance.pdist(origo))
        m_idx = np.array(np.where((dmat<max_r) & (dmat>0.01))).T
        
        # calculate scaled distances
        diff = [origo[y]-origo[x] for x,y in m_idx]
        dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
        return dotp,m_idx

    def get_other_mat(self,cutoff,atomtype):

        ii = range(len(self.sequence))
        lcs, origo = self.get_lcs(ii)
        p_pos  = []
        for j in ii:
            idx =  self.residues[j].get_idx(atomtype)
            if(idx==None):
                p_pos.append([float('nan'),float('nan'),float('nan')])
            else:
                p_pos.append(self.coords[idx])
        p_pos = np.array(p_pos)

        # prune search first                
        max_r  = np.max(definitions.f_factors)*cutoff
        dmat = distance.cdist(origo,p_pos)
        m_idx = np.array(np.where((dmat<max_r) & (dmat>0.01))).T

        # calculate scaled distances
        diff = [p_pos[y]-origo[x] for x,y in m_idx]
        dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])
        dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
        dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
        
        dotp[dotp_norm>cutoff] = 0.0

        mat = np.zeros((len(ii),len(ii),3))        
        mat[m_idx[:,0],m_idx[:,1]] = dotp
        return mat

    def get_3dmat_square(self,cutoff,ii=[]):

        ll = len(ii)
        if(ll==0):
            ll = len(self.sequence)
            ii = range(ll)            
            
        dotp,m_idx = self.get_3dmat(cutoff,ii)
        dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]
        dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
        dotp[dotp_norm>cutoff] = 0.0

        mat = np.zeros((ll,ll,3))        
        mat[m_idx[:,0],m_idx[:,1]] = dotp
        # the matrix is not rescaled!
        return mat
    

    def get_gmat(self,cutoff,ii=[]):

        ll = len(ii)

        if(ll==0):
            ll = len(self.sequence)
            ii = range(ll)            
        mat = np.zeros((ll,ll,4))
            
        dotp,m_idx = self.get_3dmat(cutoff,ii)
        # if all zeros

        if(len(dotp)==0):
            return mat

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

    def get_mat_score(self,cutoff,ii=[]):

        ll = len(ii)
        if(ll==0):
            ll = len(self.sequence)
            ii = range(ll)            
        
        dotp,m_idx = self.get_3dmat(cutoff,ii)
        dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]

        # set to zero when norm is larger than cutoff
        dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
        dotp[dotp_norm>cutoff] = 0.0
        nonzero = dotp[~np.all(dotp == 0, axis=1)]
        
        return nonzero.T

    
    # this shold not be here - but it is
    def get_annotation(self):


        lcs,origo = self.get_lcs()
        ll = len(origo)

        # hardcoded cutoff!
        cutoff = 1.58
        
        # prune search first
        max_r  = np.max(definitions.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
        pairs = []

        # dot-bracked annotation
        openings = []
        closings = []

        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            diff = origo[jj]-origo[ii]

            R1 = np.dot(diff,lcs[ii])
            R2 = np.dot(-diff,lcs[jj])

            R1_scaled = R1*definitions.scale
            R2_scaled = R2*definitions.scale
            
            D1 = np.sqrt(np.sum(R1_scaled**2))
            D2 = np.sqrt(np.sum(R2_scaled**2))
            if(D1 < cutoff and D2 < cutoff):
                
                # give generic name
                int_type = ""

                rho1 = np.sqrt(R1[0]**2 + R1[1]**2)
                rho2 = np.sqrt(R2[0]**2 + R2[1]**2)
                z1_abs = np.abs(R1[2]) 
                z2_abs = np.abs(R2[2]) 

                # check stacking 
                if(z1_abs > 2.0 and z2_abs > 2.0):

                    # rho cutoff
                    if(rho1 < 5.0 and rho2 < 5.0):
                        if(R1[2] > 2.0 and R2[2] > 2.0):
                            int_type = "><"
                            pairs.append([ii,jj,int_type])
                            continue
                        if(R1[2] > 2.0 and R2[2] < -2.0):
                            int_type = ">>"
                            pairs.append([ii,jj,int_type])
                            continue
                        if(R1[2] < -2.0 and R2[2] > 2.0):
                            int_type = "<<"
                            pairs.append([ii,jj,int_type])
                            continue
                        if(R1[2] < -2.0 and R2[2] < -2.0):
                            int_type = "<>"
                            pairs.append([ii,jj,int_type])
                            continue
                        
                        
                # check pairing
                if(z1_abs < 2.0 and z2_abs < 2.0):

                    # check  wc-pair
                    r1 = self.residues[ii].res_mytype
                    r2 = self.residues[jj].res_mytype
                    if((r1 == 'rA' and r2 == 'rU') or  (r1 == 'rU' and r2 == 'rA') or \
                           (r1 == 'rC' and r2 == 'rG') or  (r1 == 'rG' and r2 == 'rC')):
                        dev1  = R1 - definitions.wc_mean
                        dev2  = R2 - definitions.wc_mean
                        maha1 = np.einsum('...k,...kl,...l->...', dev1, definitions.inv_sigma, dev1)
                        maha2 = np.einsum('...k,...kl,...l->...', dev2, definitions.inv_sigma, dev2)
                        val1 = (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha1)
                        val2 = (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha2)
                        if(val1*val2>1.0e-08):
                            int_type = "WC"

                    # GU 
                    if((r1 == 'rU' and r2 == 'rG') or  (r1 == 'rU' and r2 == 'rG')):
                        angle1 = np.arctan2(R1[1],R1[0])
                        angle2 = np.arctan2(R2[1],R2[0])
                        if((angle1> definitions.theta1 and angle1 <= definitions.theta2) and \
                               (angle2> definitions.theta1 and angle2 <= definitions.theta2)):
                            int_type = "WC"

                    if(int_type == "WC"):
                        pairs.append([ii,jj,int_type])
                        assert (ii not in openings)
                        assert (jj not in openings)
                        assert (ii not in closings)
                        assert (jj not in closings)
                        openings.append(ii)
                        closings.append(jj)
                        continue

                # check non-canonical
                if(z1_abs < 2.0 or z2_abs < 2.0):
                    angles = [np.arctan2(R1[1],R1[0]),np.arctan2(R2[1],R2[0])]
                    for angle in angles:
                        if(angle > definitions.theta1 and angle <= definitions.theta2):
                            int_type += "W"
                        if(angle > definitions.theta2 or angle <= definitions.theta3):
                            int_type += "H"
                        if(angle <= definitions.theta1 and angle > definitions.theta3):
                            int_type += "S"
                                
                                
                if(int_type ==""):
                    int_type = "XX"
                pairs.append([ii,jj,int_type])
                # do dot-bracket annotation
                #if(int_type == "WC"):

                            
        return pairs, openings,closings

    def string_pdb(self,idx=None,noP=False):

        if(idx==None):
            idx = range(len(self.residues))

        string = "MODEL \n"
        string += "REMARK "
        string += " ".join([self.sequence_id[i] for i in idx])
        string += "\n"
        for serial,i in enumerate(idx):
            coords = self.coords[self.residues[i].first:self.residues[i].last+1]
            
            if(noP==True and serial==0):
                string += self.residues[i].pdb_string(coords,noP=True)
            else:
                string += self.residues[i].pdb_string(coords)

        string += "ENDMDL\n"
        return string



    def get_lcs_com(self,idx=[None]):

        if(idx[0]==[None]):
            idx = xrange(len(self.residues))
            
        coords1 = self.coords[self.lcs_idx[0][idx]]
        coords2 = self.coords[self.lcs_idx[1][idx]]
        coords3 = self.coords[self.lcs_idx[2][idx]]
        coords = np.array([coords1,coords2,coords3])

        # calculate center of mass
        origo = np.sum(np.sum(coords,axis=0),axis=0)/(len(idx)*3)
        return origo


