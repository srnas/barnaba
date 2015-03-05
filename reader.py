import sys
import numpy as np
from scipy.spatial import distance
#from multiprocessing import Pool,cpu_count
#from joblib import Parallel, delayed
#import profile

class Names:

    # list of all RNA atoms
    rna_atoms = ["P","OP1","OP2","O1P","O2P",\
                 "O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                 "N9","C8","N7","C5","C6","N6","N1","C2","N3","C4",\
                 "N1","C2","O2","N3","N4",\
                 "N2","O4","O6" ]

    # list of all RNA atoms
    base_atoms = ["C2","C4","C6"]
    
    # list of all RNA atoms divided per residue type. phosphate group atoms are missing (on purpose)
    atoms_l= [["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                   "N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"],\
                  ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                       "N1","C2","O2","N3","C4","N4","C5","C6"],\
                  ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                       "N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"],\
                  ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                       "N1","C2","O2","N3","C4","O4","C5","C6"]]

    rnas = ['U','rU','RU','RU5','RU3',\
                'C','rC','RC','RC5','RC3',\
                'G','rG','RG','RG5','RG3',\
                'A','rA','RA','RA5','RA3']
    rna_l = ["A","C","U","G"]
    f_factors = [5.,5.,3.]
    scale = [1./f_factors[0],1./f_factors[1],1./f_factors[2]]
    known_abbrev = ["A","C","G","U","N","Y","R","%"]

    interactions = ['..','>>','<<','<>','><','WC','WW','WS','WH','HH','HS','HW','SS','SH','SW','XX']
    pairings = ['WC','WW','WS','WH','HH','HS','HW','SS','SH','SW']
    
    # mean values and covariance matrix for wc-pair calculation
    # extracted from empirical distribution

    wc_mean = [2.86,4.67,0.01]
    wc_sigma = [[ 0.26, -0.14,0.0],\
                 [-0.14,0.13,0.0],\
                 [ 0.0, 0.0,  0.33 ]]

    inv_sigma = np.linalg.pinv(wc_sigma)
    det_sigma = np.linalg.det(wc_sigma)
    # treshold values for base pair edges were obtained
    # from the angular distribution
    
    theta1 = 0.16
    theta2 = 2.0
    theta3 = -2.0



class Atom:


    def __init__(self,atom_data):

        self.atom_n = int(atom_data[0])
        self.atom_t = atom_data[1]
        self.residue_t = atom_data[2]
        self.chain = atom_data[3]
        self.residue_n = int(atom_data[4])
        self.coords = atom_data[5:8]
    
    def pdb_string(self):
        x = self.coords[0]
        y = self.coords[1]
        z = self.coords[2]
        occ = 1.0
        bfac = 0.0
        if(len(self.atom_t)==4):
            string = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                % ("ATOM",self.atom_n,self.atom_t," ",self.residue_t,self.chain,self.residue_n," ",x,y,z,occ,bfac)
        else:
            string = "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                % ("ATOM",self.atom_n,self.atom_t," ",self.residue_t,self.chain,self.residue_n," ",x,y,z,occ,bfac)            
        return string

class Residue:

    
    def __init__(self,atoms_data):
        self.atoms = []
        self.residue_t = atoms_data[0][2]
        self.residue_n = atoms_data[0][4]
        self.chain = atoms_data[0][3]
        for i in xrange(len(atoms_data)):
            self.atoms.append(Atom(atoms_data[i]))

    def get_atom(self,atom_t):
        for i in xrange(len(self.atoms)):
            if(atom_t == self.atoms[i].atom_t):
                return self.atoms[i].coords
        return None
        
    def set_name(self,name):
        self.residue_t = name
        self.residue_n = name
        self.chain = name
        
    def pdb_string(self):
        string = ""
        for i in xrange(len(self.atoms)):
            string +=  self.atoms[i].pdb_string() 
        return string

    def get_com(self):
        center = np.array([0.0,0.0,0.0])
        for i in xrange(len(self.atoms)):
            center += self.atoms[i].coords
        cc = center/len(self.atoms)
        return cc

class Structure:

    def __init__(self,residue_data,natoms):

        self.residues = []
        self.sequence = []
        self.sequence_id = []
        self.chain = []
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
            idx = range(len(self.residues))

        coords = []
        for i in idx:
            coords.append(self.residues[i].get_com())
        return np.array(coords)

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
        for j in xrange(0,len(self.sequence)-l1):
            sub_seq = self.sequence[j:j+l1]
            sub_chain = self.chain[j:j+l1]
            if(sub_chain[0] == sub_chain[-1]):
                if(is_match("".join(sub_seq),query)): idx.append(range(j,j+l1))

        # one bulged base
        if(bulges==1):
            for j in xrange(0,len(self.sequence)-l1-1):
                sub_seq = self.sequence[j:j+l1+1]
                sub_chain = self.chain[j:j+l1+1]
                if(sub_chain[0] == sub_chain[-1]):
                    for k in xrange(1,l1):
                        tmp_idx = range(j,j+k) + range(j+k+1,j+l1+1)
                        sub_seq_tmp = [self.sequence[l] for l in tmp_idx]
                        if(is_match("".join(sub_seq_tmp),query)):idx.append(tmp_idx)

        # Two bulges
        if(bulges==2):
            for j in xrange(0,len(self.sequence)-l1-2):
                sub_seq = self.sequence[j:j+l1+2]
                sub_chain = self.chain[j:j+l1+2]
                if(sub_chain[0] == sub_chain[-1]):
                    for k in xrange(1,l1):
                        tmp_idx = range(j,j+k) + range(j+k+2,j+l1+2)
                        sub_seq_tmp = [self.sequence[l] for l in tmp_idx]
                        if(is_match("".join(sub_seq_tmp),query)):idx.append(tmp_idx)

        return idx


class Pdb:

    def __init__(self,filename,base_only=False):


        if(base_only==False):
            ok_atoms = Names.rna_atoms
        else:
            ok_atoms = Names.base_atoms

        self.models = []
        
        data = []
        natoms = []
        RES_ID_PREV = "XXX"
        num = 0

        # read file
        fh = open(filename,'r')
        for line in fh:           
            if(line[0:6].strip()=="ATOM"):

                REST = line[17:20].strip()
                RESN = line[22:26].strip()
                CHAIN = line[21:22].strip()
                ATMT = line[12:16].strip()
                ATMT.replace("*","'")

                RES_ID = REST + "_" + RESN + "_" + CHAIN

                # skip non-RNA
                if(REST not in Names.rnas):
                    continue
                # change name
                if(len(REST) > 1):
                    REST = REST[1]
                
                # skip alternate locations
                ALT = line[16:17].strip()
                if(ALT!="" and ALT !="A"):
                    err =  "# Warning: the PDB contains alternative locations for residue %s %s %s \n" % (REST,str(RESN),ALT)
                    sys.stderr.write(err)
                    continue

                # skip hydrogens and unknown atomtypes
                if(ATMT not in ok_atoms):
                    continue

                X = float(line[30:38])
                Y = float(line[38:46])
                Z = float(line[46:54])
                ATMN = line[6:11].strip()

                data.append([ATMN,ATMT,REST,CHAIN,RESN,X,Y,Z])
                if(RES_ID_PREV != RES_ID):
                    if(RES_ID_PREV!="XXX"): natoms.append(num)
                    RES_ID_PREV = RES_ID
                    num = 0
                num += 1

            if(line[0:6].strip()=="ENDMDL"):
                if(len(data) > 1):
                    natoms.append(num)
                    self.models.append(Structure(data,natoms))
                    data = []
                    natoms = []
                    num = 0
                    RES_ID_PREV = "XXX"



        fh.close()                    

        # one last time (when ENDMDL misses)
        if(len(data) > 1):
            natoms.append(num)
            self.models.append(Structure(data,natoms))



