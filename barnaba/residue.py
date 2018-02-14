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

# Make sure that range returns an iterator also in python2 (using future module)
from builtins import range

import sys
from . import definitions


class Residue:

    def __init__(self,atoms_data,idx_last):

        self.res_type = atoms_data[0][2]
        self.res_num = int(atoms_data[0][4])
        self.chain = atoms_data[0][3]
        self.res_mytype = atoms_data[0][8]
        self.res_id = atoms_data[0][9]
        self.atom_types = []
        self.atom_numbers = []
        self.atom_coords = []
        for i in range(len(atoms_data)):
            # make sure that atoms is not duplicated
            # this can happen with multiple occupancy
            if(atoms_data[i][1] in self.atom_types):
                err = "# Warning: atom %s is already present in residue %s \n" % (atoms_data[i][1],self.res_id)
                sys.stderr.write(err)
                continue
                
            self.atom_numbers.append(int(atoms_data[i][0]))
            self.atom_types.append(atoms_data[i][1])
            self.atom_coords.append([atoms_data[i][5],atoms_data[i][6],atoms_data[i][7]])
            
        self.last = idx_last
        self.first = idx_last - len(self.atom_numbers) + 1

        
    # return coordinates of given atom (accession by name)
    def __getitem__(self,atom_type):
        try:
            idx = self.atom_types.index(atom_type)
            return self.atom_coords[idx]
        except:
            err = "# Warning: no %s atom in residue %s \n" % (atom_type,self.res_id)
            sys.stderr.write(err)
            return [float('nan'),float('nan',),float('nan')]


        

        
    # return all coordinates (silent version, used in ENM)
    def get_idx(self,atom_type):
        try:
            return self.atom_types.index(atom_type) + self.first
        except:
            return float('NaN')
        


    # return all coordinates
    def get_coords(self):
        return self.atom_coords

            

    # return string with pdb
    def pdb_string(self,coords,noP=False):
        string = ""
        occ = 1.0
        bfac = 1.0

        for i in range(len(self.atom_types)):

            atype = self.atom_types[i]
            if(noP==True and atype in definitions.term_5prime):
                continue
            
            anum = self.atom_numbers[i]
            x = coords[i][0]
            y = coords[i][1]
            z = coords[i][2]
            insert = (self.res_id).split("_")[-1]
            if(insert=="0"):
                insert = " "
            if(len(atype)==4):
                string += "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                         % ("ATOM",anum,atype," ",self.res_type,self.chain,self.res_num,insert,x,y,z,occ,bfac)
            else:
                string += "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                          % ("ATOM",anum,atype," ",self.res_type,self.chain,self.res_num,insert,x,y,z,occ,bfac)
        return string
