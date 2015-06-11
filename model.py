

import numpy as np
from scipy.spatial import distance

import residue as rd
import definitions

class Model:

    def __init__(self,residues_data):

        self.sequence = []
        self.sequence_id = []
        self.residues = []
        self.chain = []
        for i in xrange(len(residues_data)):
            self.residues.append(rd.Residue(residues_data[i]))
            self.sequence.append(residues_data[i][0][2])
            self.sequence_id.append(residues_data[i][0][9])
            self.chain.append(residues_data[i][0][3])
            
    def __getitem__(self, ind):
        return self.residues[ind]

    def get_lcs(self,idxs=[None]):

        coords = []
        if(idxs[0]==None): idxs=range(len(self.sequence))
        for i in idxs:
            if(len(self.residues[i].lcs_atoms)==0):
                print self.sequence_id[i]
            else:
                coords.append(self.residues[i].lcs_atoms)

        if(len(coords)==0): return [],[]

        coords = np.array(coords)
        # calculate center of mass
        origo = np.sum(coords,axis=1)/3.0
        # CoM-C2 (x axis)
        x = coords[:,0]-origo
        x_norm = 1./np.sqrt(np.sum(x*x,axis=1))
        x = x*np.transpose([x_norm,x_norm,x_norm])
        # CoM-C4/C6 
        c = coords[:,1]-origo
        # z/y axis
        z = np.cross(x,c,axis=1)
        z_norm = 1./np.sqrt(np.sum(z*z,axis=1))
        z = z*np.transpose([z_norm,z_norm,z_norm])
        y = np.cross(z,x,axis=1)
        lcs = np.array([x.T,y.T,z.T]).T
        return lcs,origo

    def get_4dmat(self,cutoff,idxs=[None]):

        lcs,origo = self.get_lcs(idxs)

        ll = len(origo)

        fact = np.pi/cutoff
        mat_gb = np.zeros((ll,ll,4))
        
        # prune search first
        max_r  = np.max(definitions.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
        
        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            diff = origo[jj]-origo[ii]

            R1_scaled = np.dot(diff,lcs[ii])*definitions.scale
            R2_scaled = np.dot(-diff,lcs[jj])*definitions.scale
            
            D1 = np.sqrt(np.sum(R1_scaled**2))
            D2 = np.sqrt(np.sum(R2_scaled**2))
            
            
            if(D1 < cutoff):
                D1_S = D1*fact
                s = np.sin(D1_S)/D1_S
                mat_gb[ii,jj,0] = s*R1_scaled[0]
                mat_gb[ii,jj,1] = s*R1_scaled[1]
                mat_gb[ii,jj,2] = s*R1_scaled[2]
                mat_gb[ii,jj,3] = (1.0+np.cos(D1_S))/fact
                
            if(D2 < cutoff):
                D2_S = (D2*fact)
                s = np.sin(D2_S)/D2_S
                idx2 = jj*ll + ii
                mat_gb[jj,ii,0] = s*R2_scaled[0]
                mat_gb[jj,ii,1] = s*R2_scaled[1]
                mat_gb[jj,ii,2] = s*R2_scaled[2]
                mat_gb[jj,ii,3] = (1.0+np.cos(D2_S))/fact
            
        return mat_gb


    def get_mat_score(self,cutoff):
        
        lcs,origo = self.get_lcs()
        ll = len(origo)
        #mat_gb = np.zeros((ll,ll,3))
        mat = []
        
        # prune search first
        max_r  = np.max(definitions.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
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
            if(D1 < cutoff): mat.append(R1)
            if(D2 < cutoff): mat.append(R2)
                
        return np.array(mat).T

    
    def get_idx(self,query,bulges=0):
    
        def is_match(s1,s2):
            for r1,r2 in zip(s1,s2):
                if(r1==r2): continue
                if(r2=="N"): continue
                if(r2=="Y" and (r1=="C" or r1=="U")): continue
                if(r2=="R" and (r1=="G" or r1=="A")): continue
                return False
            return True

        idx = []
        l1 = len(query)
        for j in xrange(0,len(self.sequence)-l1+1):
            sub_seq = self.sequence[j:j+l1]
            sub_chain = self.chain[j:j+l1]
            if(sub_chain[0] == sub_chain[-1]):
                numbers = np.array([int(x.split("_")[0]) for x in (self.sequence_id[j:j+l1])])
                diff = numbers[1:]-numbers[:-1]
                if(any(x!=1 for x in diff) == False):
                    if(is_match("".join(sub_seq),query)):
                        idx.append(range(j,j+l1))

        # one bulged base
        if(bulges>0):
            for j in xrange(0,len(self.sequence)-l1-1):
                sub_seq = self.sequence[j:j+l1+1]
                sub_chain = self.chain[j:j+l1+1]
                if(sub_chain[0] == sub_chain[-1]):
                    for k in xrange(1,l1):
                        tmp_idx = range(j,j+k) + range(j+k+1,j+l1+1)
                        sub_seq_tmp = [self.sequence[l] for l in tmp_idx]
                        if(is_match("".join(sub_seq_tmp),query)):idx.append(tmp_idx)

        # Two bulges
        if(bulges>1):
            for j in xrange(0,len(self.sequence)-l1-2):
                sub_seq = self.sequence[j:j+l1+2]
                sub_chain = self.chain[j:j+l1+2]
                if(sub_chain[0] == sub_chain[-1]):
                    for k in xrange(1,l1):
                        tmp_idx = range(j,j+k) + range(j+k+2,j+l1+2)
                        sub_seq_tmp = [self.sequence[l] for l in tmp_idx]
                        if(is_match("".join(sub_seq_tmp),query)):idx.append(tmp_idx)

        return idx

    
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
                    r1 = self.sequence[ii]
                    r2 = self.sequence[jj]
                    if((r1 == 'A' and r2 == 'U') or  (r1 == 'U' and r2 == 'A') or \
                           (r1 == 'C' and r2 == 'G') or  (r1 == 'G' and r2 == 'C')):
                        dev1  = R1 - definitions.wc_mean
                        dev2  = R2 - definitions.wc_mean
                        maha1 = np.einsum('...k,...kl,...l->...', dev1, definitions.inv_sigma, dev1)
                        maha2 = np.einsum('...k,...kl,...l->...', dev2, definitions.inv_sigma, dev2)
                        val1 = (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha1)
                        val2 = (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha2)
                        if(val1*val2>1.0e-08):
                            int_type = "WC"

                    # GU 
                    if((r1 == 'U' and r2 == 'G') or  (r1 == 'U' and r2 == 'G')):
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

    def string_pdb(self,idx=None):

        if(idx==None):
            idx = range(len(self.residues))

        string = "MODEL \n"
        string += "REMARK "
        string += " ".join([self.sequence_id[i] for i in idx])
        string += "\n"
        for i in idx:
            string += self.residues[i].pdb_string()
        string += "ENDMDL\n"
        return string



    def get_lcs_com(self,idx=[None]):

        if(idx[0]==[None]):
            idx = xrange(len(self.residues))

        com = []
        for i in idx:
            com.extend(self.residues[i].lcs_atoms)
        return np.average(np.array(com),axis=0)


    def get_3dmat(self,cutoff):

        lcs,origo = self.get_lcs()
        ll = len(origo)
        mat_gb = np.zeros((ll,ll,3))
        
        # prune search first
        max_r  = np.max(definitions.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
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
            if(D1 < cutoff): mat_gb[ii,jj] = R1
            if(D2 < cutoff): mat_gb[jj,ii] = R2
            
        return mat_gb
