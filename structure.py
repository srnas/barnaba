class Structure:

    def __init__(self,residue_data,natoms,name):

        self.residues = []
        self.sequence = []
        self.sequence_id = []
        self.chain = []
        self.name = name
        start = 0
        for i in xrange(len(natoms)):

            res_id = residue_data[start][4] + "_"  + residue_data[start][2] + "_" + residue_data[start][3]  
            self.residues.append(Residue(residue_data[start:start+natoms[i]]))
            self.sequence.append(residue_data[start][2])

            self.sequence_id.append(res_id)
            self.chain.append(residue_data[start][3])
            start += natoms[i]
            
    def get_com(self,idx=None):

        if(idx==None):
            idx = xrange(len(self.residues))

        coords = []
        for i in idx:
            coords.append(self.residues[i].get_com())
        return np.array(coords)


    def calc_dihe(self,vecs):

        if(any(x == None for x in vecs)):
            return float('nan')
        vecs = np.array(vecs)
        b = vecs[:-1] - vecs[1:]
        norms = [np.linalg.norm(b[0]),np.linalg.norm(b[1]),np.linalg.norm(b[2])]
        if(any(x>1.7 for x in norms)):
            err = "# Warning: a bond lenght is suspiciously large: %s \n" %(self.name)
            err += "# %5.2f %5.2f %5.2f \n" % (norms[0],norms[1],norms[2])
            err += "# This dihedral will not be calculated \n"
            sys.stderr.write(err)
            return float('nan')
            
        # "Flip" the first vector so that eclipsing vectors have dihedral=0
        b[0] *= -1
        # Use dot product to find the components of b1 and b3 that are not
        # perpendicular to b2. Subtract those components. The resulting vectors
        # lie in parallel planes.
        v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
        v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
        b1 = b[1] / norms[1]
        x = np.dot(v[0], v[1])
        m = np.cross(v[0], b1)
        y = np.dot(m, v[1])

        return np.arctan2( y, x )

    def get_bb_torsions(self):
        
        idx = xrange(len(self.residues))
        angles = []
        for i in idx:
            residue = self.residues[i]

            ang = [float('nan')]*7

            at1 = residue.get_atom("P")
            at2 = residue.get_atom("O5'")
            at3 = residue.get_atom("C5'")
            at4 = residue.get_atom("C4'")
            at5 = residue.get_atom("C3'")
            at6 = residue.get_atom("O3'")            

            if(i!=0 and (residue.chain == self.residues[i-1].chain)):
                at0 = self.residues[i-1].get_atom("O3'")
                ang[0] = self.calc_dihe([at0,at1,at2,at3]) # alpha  [O3-P-O5-C5]
            ang[1] = self.calc_dihe([at1,at2,at3,at4])     # beta   [ P-O5p-C5p-C4p]
            ang[2] = self.calc_dihe([at2,at3,at4,at5])     # gamma [ O5p-C5p-C4p-C3p]
            ang[3] = self.calc_dihe([at3,at4,at5,at6])     # delta [ C5p-C4p-C3p-O3p]
                
            if(i!=idx[-1] and (residue.chain == self.residues[i+1].chain)):
                at7 = self.residues[i+1].get_atom("P")     
                at8 = self.residues[i+1].get_atom("O5'")
                ang[4] = self.calc_dihe([at4,at5,at6,at7])    # epsilon [C4p-C3p-O3p-P]
                ang[5] = self.calc_dihe([at5,at6,at7,at8])    # zeta[C3p-O3p-P-O5p]

            at9 = residue.get_atom("O4'")     
            at10 = residue.get_atom("C1'")     
            if(residue.residue_t == "A" or residue.residue_t == "G"):
                # chi (pur) [O4p-C1p-N9-C4]
                at11 = residue.get_atom("N9")
                at12 = residue.get_atom("C4")

            if(residue.residue_t == "C" or residue.residue_t == "U"):
                # chi (pyr) [O4p-C1p-N1-C2]
                at11 = residue.get_atom("N1") 
                at12 = residue.get_atom("C2")
            ang[6] = self.calc_dihe([at9,at10,at11,at12])
            angles.append(ang)

        return np.degrees(angles)


    def get_pucker(self):
        
        idx = xrange(len(self.residues))
        angles = []
        for i in idx:
            residue = self.residues[i]

            ang = [float('nan')]*7

            at1 = residue.get_atom("C4'")
            at2 = residue.get_atom("O4'")
            at3 = residue.get_atom("C1'")
            at4 = residue.get_atom("C2'")
            at5 = residue.get_atom("C3'")

            ang[0] = self.calc_dihe([at1,at2,at3,at4])      # v0: C4'-O4'-C1'-C2'
            ang[1] = self.calc_dihe([at2,at3,at4,at5])      # v1: O4'-C1'-C2'-C3'
            ang[2] = self.calc_dihe([at3,at4,at5,at1])      # v2: C1'-C2'-C3'-C4'
            ang[3] = self.calc_dihe([at4,at5,at1,at2])      # v3: C2'-C3'-C4'-O4'
            ang[4] = self.calc_dihe([at5,at1,at2,at3])      # v4: C3'-C4'-O4'-C1'


            x1=ang[4]+ang[1]-ang[3]-ang[0]
            x2=Names.Pconst*ang[2]
            p0 = np.arctan2(x1,x2) 
            if(p0<0):
                p0 += 2.0*np.pi
            tm = ang[2]/np.cos(p0)  
            ang[5] = tm
            ang[6] = p0

            angles.append(ang)

        return np.degrees(angles)
            

            
    def get_lcs(self,idx=None,permissive=True):

        if(idx==None):
            idx = range(len(self.residues))

        coords = []
        for i in idx:
            residue = self.residues[i]
            v = [residue.get_atom("C2"),residue.get_atom("C4"),residue.get_atom("C6")]
            # invert for purines!
            if(residue.residue_t == "A" or residue.residue_t == "G"):
                v = [residue.get_atom("C2"),residue.get_atom("C6"),residue.get_atom("C4")]

            # check that atom is there
            if(None in v):
                null = [0.0,0.0,0.0]
                v = [null,null,null]
                err = "# Warning: No C2,C4 or C6 atoms in residue %s %s \n" % (residue.residue_n, residue.residue_t)
                sys.stderr.write(err)
                if(permissive==False):
                    return None,None
                
            coords.append(v)
        if(len(coords)<2):
            err =  "# Warning: less than 2 RNA residues in PDB \n"
            sys.stderr.write(err)
            return None,None
        
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



    def get_4dmat(self,cutoff,idx=None,permissive=True):


        lcs,origo = self.get_lcs(idx,permissive)
        if(lcs == None and permissive==False):
            return None

        #print lcs
        ll = len(origo)

        fact = np.pi/cutoff
        mat_gb = np.zeros((ll,ll,4))
        
        # prune search first
        max_r  = np.max(Names.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
        
        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            diff = origo[jj]-origo[ii]

            R1_scaled = np.dot(diff,lcs[ii])*Names.scale
            R2_scaled = np.dot(-diff,lcs[jj])*Names.scale
            
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



    def get_3dmat(self,cutoff):

        lcs,origo = self.get_lcs()
        ll = len(origo)
        mat_gb = np.zeros((ll,ll,3))
        
        # prune search first
        max_r  = np.max(Names.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            diff = origo[jj]-origo[ii]

            R1 = np.dot(diff,lcs[ii])
            R2 = np.dot(-diff,lcs[jj])

            R1_scaled = R1*Names.scale
            R2_scaled = R2*Names.scale
            
            D1 = np.sqrt(np.sum(R1_scaled**2))
            D2 = np.sqrt(np.sum(R2_scaled**2))
            if(D1 < cutoff): mat_gb[ii,jj] = R1
            if(D2 < cutoff): mat_gb[jj,ii] = R2
            
        return mat_gb

    def get_mat_score(self,cutoff):

        lcs,origo = self.get_lcs()
        ll = len(origo)
        #mat_gb = np.zeros((ll,ll,3))
        mat = []
        
        # prune search first
        max_r  = np.max(Names.f_factors)*cutoff
        dmat = distance.pdist(origo)
        c_idx = (dmat<max_r).nonzero()[0]
        m_idx = np.array(np.triu_indices(ll,1)).T[c_idx] 
        
        for kk in xrange(len(c_idx)):
            ii = m_idx[kk][0]
            jj = m_idx[kk][1]
            diff = origo[jj]-origo[ii]

            R1 = np.dot(diff,lcs[ii])
            R2 = np.dot(-diff,lcs[jj])

            R1_scaled = R1*Names.scale
            R2_scaled = R2*Names.scale
            
            D1 = np.sqrt(np.sum(R1_scaled**2))
            D2 = np.sqrt(np.sum(R2_scaled**2))
            if(D1 < cutoff): mat.append(R1)
            if(D2 < cutoff): mat.append(R2)

        return np.array(mat).T


        

    def get_annotation(self):


        lcs,origo = self.get_lcs()
        ll = len(origo)

        # hardcoded cutoff!
        cutoff = 1.58
        
        # prune search first
        max_r  = np.max(Names.f_factors)*cutoff
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

            R1_scaled = R1*Names.scale
            R2_scaled = R2*Names.scale
            
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
                        dev1  = R1 - Names.wc_mean
                        dev2  = R2 - Names.wc_mean
                        maha1 = np.einsum('...k,...kl,...l->...', dev1, Names.inv_sigma, dev1)
                        maha2 = np.einsum('...k,...kl,...l->...', dev2, Names.inv_sigma, dev2)
                        val1 = (2 * np.pi)**(-1.5)* Names.det_sigma**(0.5)*np.exp(-0.5 * maha1)
                        val2 = (2 * np.pi)**(-1.5)* Names.det_sigma**(0.5)*np.exp(-0.5 * maha2)
                        if(val1*val2>1.0e-08):
                            int_type = "WC"

                    # GU 
                    if((r1 == 'U' and r2 == 'G') or  (r1 == 'U' and r2 == 'G')):
                        angle1 = np.arctan2(R1[1],R1[0])
                        angle2 = np.arctan2(R2[1],R2[0])
                        if((angle1> Names.theta1 and angle1 <= Names.theta2) and \
                               (angle2> Names.theta1 and angle2 <= Names.theta2)):
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
                        if(angle > Names.theta1 and angle <= Names.theta2):
                            int_type += "W"
                        if(angle > Names.theta2 or angle <= Names.theta3):
                            int_type += "H"
                        if(angle <= Names.theta1 and angle > Names.theta3):
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
