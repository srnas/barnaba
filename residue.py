import sys

pyr = ["rC","rU","dC","dU","dT"]
pur = ["rA","rG","dA","dG"]

class Residue:

    def __init__(self,atoms_data):

        self.res_type = atoms_data[0][2]
        self.res_num = int(atoms_data[0][4])
        self.chain = atoms_data[0][3]
        #self.mol_id = atoms_data[0][8]
        self.res_mytype = atoms_data[0][8]
        self.res_id = atoms_data[0][9]
        self.atom_types = []
        self.atom_coords = []
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
            self.atom_coords.append([atoms_data[i][5],atoms_data[i][6],atoms_data[i][7]])
        #self.lcs_atoms = self.get_lcs_coords()

        
    # return coordinates of given atom (accession by name)
    def __getitem__(self,atom_type):
        try:
            idx = self.atom_types.index(atom_type)
            return self.atom_coords[idx]
        except:
            err = "# Warning: no %s atom in residue %s \n" % (atom_type,self.res_id)
            sys.stderr.write(err)
            return [-9999.,-9999.,-9999.]


    # return all coordinates (silent version, used in ENM)
    def get_atom(self,atom_type):
        try:
            idx = self.atom_types.index(atom_type)
            return self.atom_coords[idx]
        except:
            return []

    # return all coordinates
    def get_coords(self):
        return self.atom_coords

    # return all lcs coordinates
    def get_lcs_coords(self):
        
        #mycoords = []
        vec = []

        if(self.res_mytype in pyr):
            vec = [self["C2"],self["C4"],self["C6"]]
        if(self.res_mytype in pur):
            vec = [self["C2"],self["C6"],self["C4"]]
        return vec
            


    # return string with pdb
    def pdb_string(self):
        string = ""
        occ = 1.0
        bfac = 1.0

        for i in xrange(len(self.atom_types)):

            atype = self.atom_types[i]
            anum = self.atom_numbers[i]
            x = self.atom_coords[i][0]
            y = self.atom_coords[i][1]
            z = self.atom_coords[i][2]
            
            if(len(atype)==4):
                string += "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                         % ("ATOM",anum,atype," ",self.res_type,self.chain,self.res_num," ",x,y,z,occ,bfac)
            else:
                string += "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
                          % ("ATOM",anum,atype," ",self.res_type,self.chain,self.res_num," ",x,y,z,occ,bfac)
        return string
