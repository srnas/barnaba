#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2017 Sandro Bottaro (sandro.bottaro@bio.ku.dk)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import, division, print_function
import sys
import numpy as np
from . import definitions
from . import model as md

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

        
    residue_dict = {'U': 'U', 'rU':'U','RU':'U','RU5':'U','RU3':'U','U3':'U','U5':'U',\
                        'H2U':'U','PSU':'U','OMU':'U','UR3':'U','5MU':'U',\
                        'A': 'A', 'rA':'A','RA':'A','RA5':'A','RA3':'A','A3':'A','A5':'A',\
                        '1MA':'A',\
                        'C': 'C', 'rC':'C','RC':'C','RC5':'C','RC3':'C','C3':'C','C5':'C',\
                        'OMC':'C','5MC':'C',\
                        'G': 'G', 'rG':'G','RG':'G','RG5':'G','RG3':'G','G3':'G','G5':'G',\
                        '2MG':'G','YG':'G','7MG':'G','OMG':'G','1MG':'G',\
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

    def __init__(self,filename,res_mode,permissive=False):

        self.ok_residues = []
        if("R" in res_mode):
            self.ok_residues.extend(Names.rna_residues[:])
        if("S" in res_mode):
            self.ok_residues.extend(Names.rna_special[:])
        if("D" in res_mode):
            self.ok_residues.extend(Names.dna_residues[:])
        if("P" in res_mode):
            self.ok_residues.extend(Names.prt_residues[:])
        self.time = 0
        self.filename = filename
        self.fh = open(filename,'r')
        self.natoms = 0
        self.model = None
        print("# Initializing file", filename)
        
        self.parse(permissive)


    def parse(self,permissive=False):
        
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
            
            # not used - but in case
            #segid = line[72:76]
            #bfactor = float(line[60:66])
            #occupancy = float(line[54:60])
            insertion = line[26]
            if(insertion==" "):
                insertion = "0"
            res_mytype = Names.residue_dict[res_type]
            res_id = "%s_%s_%s_%s" % (res_num,res_type,chain,insertion)
            
            vv = [atom_num,atom_type,res_type,chain,res_num,X,Y,Z,res_mytype,res_id]
            return vv

        # read file
        tmp_data = []
        for line in self.fh:
            at = line[0:6].strip()
            if(at =="ATOM" or at == "HETATM"):
                # skip residues
                res_type = line[17:20].strip()
                if(res_type not in self.ok_residues): continue
                
                # skip if alternative position is different from A or nothing
                if(line[16]  != " " and line[16] != "A"): 
                    sys.stderr.write("# Warning: skipping residue with multiple occupancy %s \n" % (line[17:26]))
                    continue
                # skip insertions as well
                #insertion = line[26]
                #if(insertion!=" "):
                #    sys.stderr.write("# Warning: skipping insertion residue %s \n" % (line[17:27]))
                #    continue

                #print(readline(line))
                tmp_data.append(readline(line))
                self.natoms += 1
            if(at=="ENDMDL"):
                self.model = md.Model(tmp_data,not permissive)
                return 0
        
        if(self.model==None):
            self.model = md.Model(tmp_data,not permissive)
            return 0
            
                

        
    def read(self):

        assert self.model != None, "# Model uninitialized. call parse before read!"

        # read pdb
        data = []
        for line in self.fh:
            at = line[0:6].strip()
            if(at =="ATOM" or at == "HETATM"):
                # skip residues
                res_type = line[17:20].strip()
                if(res_type not in self.ok_residues): continue
                
                # skip if alternative position is different from A or nothing
                if(line[16]  != " " and line[16] != "A"): 
                    sys.stderr.write("# Warning: skipping residue with multiple occupancy %s \n" % (line[17:26]))
                    continue

                vv = [float(line[30:38]),float(line[38:46]) ,float(line[46:54])]
                data.append(vv)
            if(at=="ENDMDL"):
                self.model.set_coords(np.array(data))
                self.time += 1
                # return step
                return self.time
        self.fh.close()
        return -1
        

