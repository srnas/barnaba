import definitions
import model as md
import sys

class Names:


    rna_residues = ["U","rU","RU","RU5","RU3","U3","U5",\
                    "C","rC","RC","RC5","RC3","C3","C5",\
                    "G","rG","RG","RG5","RG3","G3","G5",\
                    "A","rA","RA","RA5","RA3","A3","A5"]
    
    rna_special = ["2MG","H2U","OMC","YG","PSU","5MC","7MG","1MA","OMU","OMG","UR3","1MG","5MU"]

    dna_residues = ["T","dT","DT","DT5","DT3","dC","DC","DC5","DC3",\
                    "dG","DG","DG5","DG3","dA","DA","DA5","DA3"]
    prt_residues = ["ALA","ARG","ASN","ASP","ASPP","CYS","GLN","GLU","GLY",\
                    "HSD","HSE","HSP","ILE","LEU","LYS","LSN","MET","PHE",\
                    "PRO","SER","THR","TRP","TYR","VAL"]

    ions = ["MG","HOH"]
    
    known_abbrev = ["A","C","G","U","N","Y","R","%"]

        
    residue_dict = {'U': 'rU', 'rU':'rU','RU':'rU','RU5':'rU','RU3':'rU','U3':'rU','U5':'rU',\
                        'H2U':'rU','PSU':'rU','OMU':'rU','UR3':'rU','5MU':'rU',\
                        'A': 'rA', 'rA':'rA','RA':'rA','RA5':'rA','RA3':'rA','A3':'rA','A5':'rA',\
                        '1MA':'rA',\
                        'C': 'rC', 'rC':'rC','RC':'rC','RC5':'rC','RC3':'rC','C3':'rC','C5':'rC',\
                        'OMC':'rC','5MC':'rC',\
                        'G': 'rG', 'rG':'rG','RG':'rG','RG5':'rG','RG3':'rG','G3':'rG','G5':'rG',\
                        '2MG':'rG','YG':'rG','7MG':'rG','OMG':'rG','1MG':'rG',\
                        'T': 'dT', 'dT':'dT','DT':'dT','DT5':'dT','DT3':'dT',\
                        'dA':'dA','DA':'dA','DA5':'dA','DA3':'dA',\
                        'dC':'dC','DC':'dC','DC5':'dC','DC3':'dC',\
                        'dG':'dG','DG':'dG','DG5':'dG','DG3':'dG',\
                        "ALA":"ALA","ARG":"ARG","ASN":"ASN","ASP":"ASP",\
                        "CYS":"CYS","GLN":"GLN","GLU":"GLU","GLY":"GLY",\
                        "HSD":"HSD","HSE":"HSD","HSP":"HSD","ILE":"ILE",\
                        "LEU":"LEU","LYS":"LYS","LSN":"LYS","MET":"MET",\
                        "PHE":"PHE","PRO":"PRO","SER":"SER","THR":"THR",\
                        "TRP":"TRP","TYR":"TYR","VAL":"VAL","HOH":"HOH"};

class Pdb:

    def __init__(self,filename,res_mode,at_mode,verbose=False):

        ok_residues = []
        if("R" in res_mode):
            ok_residues.extend(Names.rna_residues[:])
        if("S" in res_mode):
            ok_residues.extend(Names.rna_special[:])
        if("D" in res_mode):
            ok_residues.extend(Names.dna_residues[:])
        if("P" in res_mode):
            ok_residues.extend(Names.prt_residues[:])

        self.models = []
        self._parse(filename,ok_residues,at_mode,verbose)

        
    def _parse(self,filename,ok_residues,at_mode,verbose):

        def readline(line):
            res_type = line[17:20].strip()
            res_num = line[22:26].strip()
            atom_type = line[12:16].strip()
            atom_num = line[6:11].strip()

            chain = line[21:22].strip()
            alt = line[16].strip()
            X = float(line[30:38])
            Y = float(line[38:46])
            Z = float(line[46:54])
            insertion = line[26]

            # not used - but in case
            segid = line[72:76]
            bfactor = float(line[60:66])
            occupancy = float(line[54:60])
            res_mytype = Names.residue_dict[res_type]
            res_id = res_num + "_" + res_type + "_" + chain + "_" + insertion
            #mol_id = "XXX"
            #if(res_type in Names.rna_residues): mol_id = "RNA"
            #if(res_type in Names.dna_residues): mol_id = "DNA"
            #if(res_type in Names.prt_residues): mol_id = "PROT"
            
            vv = [atom_num,atom_type,res_type,chain,res_num,X,Y,Z,res_mytype,res_id]
            return vv
        
        # read file
        models_tmp  = []
        model = []
        print "# Reading", filename
        fh = open(filename,'r')
        for line in fh:
            at = line[0:6].strip()
            
            if(at=="ENDMDL" and len(model)!=0):
                models_tmp.append(model)
                model = []
                
            if(at =="ATOM" or at == "HETATM"):

                # skip residues
                res_type = line[17:20].strip()
                #print res_type
                if(res_type not in ok_residues): continue
                #print res_type

                # skip atoms
                atom_type = line[12:16].strip()
                if(at_mode == "ALL"  and "H" in atom_type): continue
                if(at_mode == "LCS" and atom_type not in definitions.rna_lcs): continue
                if(at_mode == "PUCKER" and atom_type not in definitions.rna_pucker): continue
                if(at_mode == "BB" and atom_type not in definitions.rna_torsion): continue

                # skip if alternative position is different from A or nothing
                if(line[16]  != " " and line[16] != "A"): 
                    if(verbose): sys.stderr.write("# Warning: skipping residue with multiple occupancy %s \n" % (line[16]))
                    continue


                model.append(readline(line))
                
        fh.close()
        # one last time if ENDMDL is missing
        if(len(model)!=0):
            models_tmp.append(model)
            
        # make sure there is at least one residue in file
        if(len(models_tmp)==0):
            err = "# Error: no valid residues in PDB file %s \n" % filename
            sys.stderr.write(err)

        # rearrange now
        self.models = []
        for i in xrange(len(models_tmp)):
            data_tmp = []
            mod = models_tmp[i]
            model = []
            for j in xrange(len(mod)):

                data_tmp.append(mod[j])

                cur_id = mod[j][-1]
                if(j+1==len(mod)): next_id = "XXX"
                else: next_id = mod[j+1][-1]

                if(next_id != cur_id):
                    model.append(data_tmp)
                    data_tmp = []
            self.models.append(md.Model(model))
        
            
                

        
        
