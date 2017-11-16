import definitions
import numpy as np
import sys

class Nucleic:

    def __init__(self,topology):

        
        # loop over residues
        ok_residues = []
        indeces_lcs = []
        indeces_glyco = []
        
        #self.bonds = list(topology.bonds)
        self.rna_seq_id = []
        self.rna_seq = []
        
        for res in topology.residues:
            if(res.name in definitions.residue_dict):
                res_type = definitions.residue_dict[res.name]

            else:
                if(res.name in definitions.modified_dict):
                    res_type = definitions.modified_dict[res.name]
                    warn = "# Treating nucleotide %s as %s \n" % (res,res_type)
                    sys.stderr.write(warn)
                else:
                    if(res.name not in definitions.others):
                        warn = "# Skipping unknown residue %s \n" % res 
                        sys.stderr.write(warn)
                    continue
            # try to fetch the fundamental atoms: C2,C4,C6,C1' and N1/N9
            nat = "N1"
            if res_type in definitions.purines:
                nat = "N9"
            try:
                i0 = res.atom("C2")
                i1 = res.atom("C6")
                i2 = res.atom("C4")
                i3 = res.atom("C1'")
                i4 = res.atom(nat)
            except:
                warn = "# Skipping residue %s -  missing atoms \n" % res 
                sys.stderr.write(warn)
                continue
        
            indeces_lcs.append([i0.index,i1.index,i2.index])
            indeces_glyco.append([i4.index,i3.index])
            ok_residues.append(res)
            self.rna_seq_id.append(res_type)
            my_resname = "%s_%s_%s" % (res.name,res.resSeq,res.chain.index)
            self.rna_seq.append(my_resname)
            
        self.indeces_lcs = np.asarray(indeces_lcs).T
        self.indeces_glyco = np.asarray(indeces_glyco).T
        self.ok_residues = ok_residues

        if(len(self.ok_residues)<1):
            warn = "# Only %d  found in structure. Exiting \n" % len(ok_residues) 
            sys.stderr.write(warn)
            sys.exit(1)
        #else:
        #    warn = "# %d nucleotides found in structure \n" % len(ok_residues) 
        #    sys.stderr.write(warn)
            


    
    def get_bb_torsion_idx(self, residues=None):

        ll = len(self.ok_residues)
        # if nothing is specified, all residues are considered
        if(residues!=None):
            idx_residues = []
            for k in range(len(residues)):
                if(residues[k] in self.rna_seq):
                    idx_residues.append((self.rna_seq).index(residues[k]))
                else:
                    msg = "# Fatal error. Residue \"%s\" not found \n" % residues[k]
                    msg += "# The list of residues is: %s \n" % self.rna_seq
                    sys.stderr.write(msg)
                    sys.exit(1)
        else:
            idx_residues = np.arange(ll)

        idxs = np.zeros((len(idx_residues),7,4),dtype=int)
        rr = []

        for j,i in enumerate(idx_residues):

            res = self.ok_residues[i]
            rr.append(self.rna_seq[i])
            # alpha
            try:
                res_minus = self.ok_residues[i-1]
                #print res,res_minus, res_minus.chain.index,res.chain.index
                idxs_tmp  = [res_minus.atom("O3'").index,res.atom("P").index,res.atom("O5'").index,res.atom("C5'").index]
                if(i!=0 and (res_minus.chain.index == res.chain.index)):
                    #print "mm", idxs_tmp
                    idxs[j,0] = idxs_tmp
                else:
                    pass
                
            except:
                pass
            # beta
            try: idxs[j,1] = [res.atom("P").index, res.atom("O5'").index, res.atom("C5'").index, res.atom("C4'").index]
            except: pass

            # gamma
            try: idxs[j,2] = [res.atom("O5'").index, res.atom("C5'").index, res.atom("C4'").index, res.atom("C3'").index]
            except: pass

            # delta
            try: idxs[j,3] = [res.atom("C5'").index, res.atom("C4'").index, res.atom("C3'").index, res.atom("O3'").index]
            except: pass

            # epsilon
            try:
                res_plus = self.ok_residues[i+1]
                idxs_tmp  = [res.atom("C4'").index,res.atom("C3'").index,res.atom("O3'").index,res_plus.atom("P").index]
                if(res_plus.chain.index == res.chain.index):
                   idxs[j,4] = idxs_tmp
                else:
                    pass
            except:
                pass
            
            # zeta
            try:
                res_plus = self.ok_residues[i+1]
                idxs_tmp  = [res.atom("C3'").index,res.atom("O3'").index,res_plus.atom("P").index,res_plus.atom("O5'").index]
                if(res_plus.chain.index == res.chain.index):
                   idxs[j,5] = idxs_tmp
                else:
                    pass
            except:
                pass

            # chi
            try:
                # try to find N9 first
                idxs[j,6] = [res.atom("O4'").index, res.atom("C1'").index, res.atom("N9").index, res.atom("C4").index]
            except:
                try: idxs[j,6] = [res.atom("O4'").index, res.atom("C1'").index, res.atom("N1").index, res.atom("C2").index]
                except: pass

        return idxs, rr

    def get_sugar_torsion_idx(self,residues=None):

        ll = len(self.ok_residues)
        
        # if nothing is specified, all residues are considered

        if(residues!=None):
            idx_residues = []
            for k in range(len(residues)):
                if(residues[k] in self.rna_seq):
                    idx_residues.append((self.rna_seq).index(residues[k]))
                else:
                    msg = "# Fatal error. Residue \"%s\" not found \n" % residues[k]
                    msg += "# The list of residues is: %s \n" % self.rna_seq
                    sys.stderr.write(msg)
                    sys.exit(1)
        else:
            idx_residues = np.arange(ll)


        idxs = np.zeros((len(idx_residues),5,4),dtype=int)
        rr = []

        for j,i in enumerate(idx_residues):

            res = self.ok_residues[i]
            rr.append(self.rna_seq[i])

            try: idxs[j,0]  = [res.atom("C4'").index, res.atom("O4'").index, res.atom("C1'").index, res.atom("C2'").index]
            except: pass
            
            try: idxs[j,1]  = [res.atom("O4'").index, res.atom("C1'").index, res.atom("C2'").index, res.atom("C3'").index]
            except: pass
            
            try: idxs[j,2]  = [res.atom("C1'").index, res.atom("C2'").index, res.atom("C3'").index, res.atom("C4'").index]
            except: pass
            
            try: idxs[j,3]  = [res.atom("C2'").index, res.atom("C3'").index, res.atom("C4'").index, res.atom("O4'").index]
            except: pass
            
            try: idxs[j,4]  = [res.atom("C3'").index, res.atom("C4'").index, res.atom("O4'").index, res.atom("C1'").index]
            except: pass
            
        return idxs, rr


    def get_coupling_idx(self,residues=None):

        ll = len(self.ok_residues)
        
        # if nothing is specified, all residues are considered

        if(residues!=None):
            idx_residues = []
            for k in range(len(residues)):
                if(residues[k] in self.rna_seq):
                    idx_residues.append((self.rna_seq).index(residues[k]))
                else:
                    msg = "# Fatal error. Residue \"%s\" not found \n" % residues[k]
                    msg += "# The list of residues is: %s \n" % self.rna_seq
                    sys.stderr.write(msg)
                    sys.exit(1)
        else:
            idx_residues = np.arange(ll)


        idxs = np.zeros((len(idx_residues),7,4),dtype=int)
        rr = []
        for j,i in enumerate(idx_residues):

            res = self.ok_residues[i]
            rr.append(self.rna_seq[i])
            # sugar couplings 

            try:
                idxs[j,0]  = [res.atom("H1'").index, res.atom("C1'").index, res.atom("C2'").index, res.atom("1H2'").index]
            except:
                try:
                    idxs[j,0]  = [res.atom("H1'").index, res.atom("C1'").index, res.atom("C2'").index, res.atom("H2'").index]
                except:
                    pass
            
            try:
                idxs[j,1]  = [res.atom("1H2'").index, res.atom("C2'").index, res.atom("C3'").index, res.atom("H3'").index]
            except:
                try:
                    idxs[j,1]  = [res.atom("H2'").index, res.atom("C2'").index, res.atom("C3'").index, res.atom("H3'").index]
                except:
                    pass
            
            try:
                idxs[j,2]  = [res.atom("H3'").index, res.atom("C3'").index, res.atom("C4'").index, res.atom("H4'").index]
            except:
                pass
            # beta
            try: idxs[j,3] = [res.atom("P").index, res.atom("O5'").index, res.atom("C5'").index, res.atom("C4'").index]
            except: pass

            # gamma
            try: idxs[j,4] = [res.atom("O5'").index, res.atom("C5'").index, res.atom("C4'").index, res.atom("C3'").index]
            except: pass

            # epsilon
            try:
                res_plus = self.ok_residues[i+1]
                idxs_tmp  = [res.atom("C4'").index,res.atom("C3'").index,res.atom("O3'").index,res_plus.atom("P").index]
                if(res_plus.chain.index == res.chain.index):
                   idxs[j,5] = idxs_tmp
                else:
                    pass
            except:
                pass
            
            # chi
            if(self.rna_seq_id[i]=="A" or self.rna_seq_id[i]=="G"):

                try: idxs[j,6] = [res.atom("O4'").index, res.atom("C1'").index, res.atom("N9").index, res.atom("C4").index]
                except: pass
            else:
                try: idxs[j,6] = [res.atom("O4'").index, res.atom("C1'").index, res.atom("N1").index, res.atom("C2").index]
                except: pass
                    
        return idxs, rr

 
