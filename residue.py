import sys
import definitions


class Residue:

    def __init__(self,atoms_data,idx_last):

        self.res_type = atoms_data[0][2]
        self.res_num = int(atoms_data[0][4])
        self.chain = atoms_data[0][3]
        self.res_mytype = atoms_data[0][8]
        self.res_id = atoms_data[0][9]
        self.atom_types = []
        self.atom_numbers = []
        for i in xrange(len(atoms_data)):
            # make sure that atoms is not duplicated
            # this can happen with multiple occupancy
            if(atoms_data[i][1] in self.atom_types):
                err = "# Warning: atom %s is already present in residue %s \n" % (atoms_data[i][1],self.res_id)
                sys.stderr.write(err)
                continue
                
            self.atom_numbers.append(int(atoms_data[i][0]))
            self.atom_types.append(atoms_data[i][1])
            
            #self.atom_coords.append([atoms_data[i][5],atoms_data[i][6],atoms_data[i][7]])
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
            return None
        

    # return all coordinates (silent version, used in ENM)
    #def get_atom(self,atom_type):
    #    try:
    #        idx = self.atom_types.index(atom_type)
    #        return self.atom_coords[idx]
    #    except:
    #        return []

    # return all coordinates
    def get_coords(self):
        return self.atom_coords

            


    # return string with pdb
    def pdb_string(self,coords,noP=False):
        string = ""
        occ = 1.0
        bfac = 1.0

        for i in xrange(len(self.atom_types)):

            atype = self.atom_types[i]
            if(noP==True and atype in definitions.term_5prime):
                continue
            
            anum = self.atom_numbers[i]
            x = coords[i][0]
            y = coords[i][1]
            z = coords[i][2]
            
            if(len(atype)==4):
                string += "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                         % ("ATOM",anum,atype," ",self.res_type,self.chain,self.res_num," ",x,y,z,occ,bfac)
            else:
                string += "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                          % ("ATOM",anum,atype," ",self.res_type,self.chain,self.res_num," ",x,y,z,occ,bfac)
        return string
