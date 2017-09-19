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
            if res_type in definitions.purines:
                try:
                    i0 = res.atom("C2")
                    i1 = res.atom("C6")
                    i2 = res.atom("C4")
                    i3 = res.atom("C1'")
                    i4 = res.atom("N9")
                except:
                    warn = "# Skipping residue %s -  missing atoms \n" % res 
                    sys.stderr.write(warn)
                    continue
            else:
                try:
                    i0 = res.atom("C2")
                    i1 = res.atom("C4")
                    i2 = res.atom("C6")
                    i3 = res.atom("C1'")
                    i4 = res.atom("N1")
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
            


    def get_couplings_torsion_idx(self):

        idxs = np.zeros((len(self.ok_residues),len(definitions.j3),4),dtype=int)
        for i,res in enumerate(self.ok_residues):
            for k in range(len(definitions.j3)):
                
                # take next P atoms for C4Pe and H3P couplings
                if(definitions.j3[k][0] == "C4Pe" or definitions.j3[k][0]=="H3P"):
                    try:
                        res_plus = self.ok_residues[i+1]
                        vv = [res.atom(definitions.j3[k][1][l]).index for l in range(3)]
                        vv.append(res_plus.atom(definitions.j3[k][1][3]).index)
                        idxs[i,k] = vv
                    except:
                        pass
                else:
                    try:
                        idxs[i,k] = [res.atom(at).index for at in definitions.j3[k][1]]
                    except:
                        pass

        return idxs
    
    def get_bb_torsion_idx(self):


        idxs = np.zeros((len(self.ok_residues),7,4),dtype=int)
        for i,res in enumerate(self.ok_residues):

            # alpha
            try:
                res_minus = self.ok_residues[i-1]
                #print res,res_minus, res_minus.chain.index,res.chain.index
                idxs_tmp  = [res_minus.atom("O3'").index,res.atom("P").index,res.atom("O5'").index,res.atom("C5'").index]
                if(i!=0 and (res_minus.chain.index == res.chain.index)):
                    #print "mm", idxs_tmp
                    idxs[i,0] = idxs_tmp
                else:
                    pass
                
            except:
                pass
            # beta
            try: idxs[i,1] = [res.atom("P").index, res.atom("O5'").index, res.atom("C5'").index, res.atom("C4'").index]
            except: pass

            # gamma
            try: idxs[i,2] = [res.atom("O5'").index, res.atom("C5'").index, res.atom("C4'").index, res.atom("C3'").index]
            except: pass

            # delta
            try: idxs[i,3] = [res.atom("C5'").index, res.atom("C4'").index, res.atom("C3'").index, res.atom("O3'").index]
            except: pass

            # epsilon
            try:
                res_plus = self.ok_residues[i+1]
                idxs_tmp  = [res.atom("C4'").index,res.atom("C3'").index,res.atom("O3'").index,res_plus.atom("P").index]
                if(res_plus.chain.index == res.chain.index):
                   idxs[i,4] = idxs_tmp
                else:
                    pass
            except:
                pass
            
            # zeta
            try:
                res_plus = self.ok_residues[i+1]
                idxs_tmp  = [res.atom("C3'").index,res.atom("O3'").index,res_plus.atom("P").index,res_plus.atom("O5'").index]
                if(res_plus.chain.index == res.chain.index):
                   idxs[i,5] = idxs_tmp
                else:
                    pass
            except:
                pass

            # chi
            try:
                # try to find N9 first
                idxs[i,6] = [res.atom("O4'").index, res.atom("C1'").index, res.atom("N9").index, res.atom("C4").index]
            except:
                try: idxs[i,6] = [res.atom("O4'").index, res.atom("C1'").index, res.atom("N1").index, res.atom("C2").index]
                except: pass
        return idxs

    def get_sugar_torsion_idx(self):
            
        idxs = np.zeros((len(self.ok_residues),5,4),dtype=int)
        for i,rr in enumerate(self.ok_residues):
            try: idxs[i,0]  = [rr.atom("C4'").index, rr.atom("O4'").index, rr.atom("C1'").index, rr.atom("C2'").index]
            except: pass
            
            try: idxs[i,1]  = [rr.atom("O4'").index, rr.atom("C1'").index, rr.atom("C2'").index, rr.atom("C3'").index]
            except: pass
            
            try: idxs[i,2]  = [rr.atom("C1'").index, rr.atom("C2'").index, rr.atom("C3'").index, rr.atom("C4'").index]
            except: pass
            
            try: idxs[i,3]  = [rr.atom("C2'").index, rr.atom("C3'").index, rr.atom("C4'").index, rr.atom("O4'").index]
            except: pass
            
            try: idxs[i,4]  = [rr.atom("C3'").index, rr.atom("C4'").index, rr.atom("O4'").index, rr.atom("C1'").index]
            except: pass
        return idxs
