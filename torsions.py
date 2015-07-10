#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014 Sandro Bottaro (sbottaro@sissa.it)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import reader as reader
import numpy as np
import definitions


def dihedral(b0,b1,b2,norm):

    ss = np.sum(b0*b1)
    v0= b0-((np.sum(b0*b1)*b1)*norm)
    v2= b2-((np.sum(b2*b1)*b1)*norm)
    x = np.sum(v0*v2)
    m = np.cross(v0,b1)*np.sqrt(norm)
    y = np.sum(m*v2)
    return np.arctan2( y, x )

    

def bb_angles(coords):

    diffs = coords[:-1]-coords[1:]

    norm_sq = np.sum(diffs**2,axis=1)
    
    # set to zero when atom is missing
    not_me = np.where(norm_sq>definitions.maxbond_sq)[0]
    no_idx = []
    for i in not_me:
        if(i>1): no_idx.append(i-1)
        no_idx.append(i)
        if(i<norm_sq.shape[0]-1): no_idx.append(i+1)
    norm_sq[no_idx] = float('nan')
    norm_sq_inv = 1./norm_sq

    angles = [float('nan')] # First is alpha 
    for i in xrange(coords.shape[0]-3):
        angles.append(dihedral(-diffs[i],diffs[i+1],diffs[i+2],norm_sq_inv[i+1]))
    
    angles.append(float('nan')) # Add two more angles at the end!
    angles.append(float('nan')) # Add two more angles at the end!

    return np.array(angles).reshape(-1,6)

def chi_angles(coords):
    
    diffs = coords[:,:-1]-coords[:,1:]

    norm_sq = np.sum(diffs**2,axis=2)
    norm_sq[np.where(norm_sq>definitions.maxbond_sq)[0],:] = 0.0

    b0 = -diffs[:,0]
    b1 = diffs[:,1]
    b2 = diffs[:,2]
    ss = np.sum(b0*b1,axis=1)

    v0= b0-((np.sum(b0*b1,axis=1)*b1.T)/norm_sq[:,1]).T
    v2= b2-((np.sum(b2*b1,axis=1)*b1.T)/norm_sq[:,1]).T
    x = np.sum(v0*v2,axis=1)
    m = np.cross(v0,b1).T/np.sqrt(norm_sq[:,1])
    y = np.sum(m.T*v2,axis=1)
    
    return  np.arctan2(y,x)



    #return np.degrees(angles.T),idx_up,idx_low
#
def torsions(args):


    print "# Calculating torsion angles..."
    files = args.files
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    for i in xrange(0,len(files)):
        
        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode,at_mode="TOR")

        for j in xrange(len(cur_pdb.models)):

            mod = cur_pdb.models[j]

            # fetch bb atoms and chi atoms
            coords_bb = []
            coords_chi = []
            for k in xrange(len(mod.residues)):

                # get atoms for backbone
                cc_tmp = [mod.residues[k][at]  for at in definitions.rna_backbone]
                coords_bb.extend(cc_tmp)

                # get atoms for chi angles
                cc_tmp = []
                atoms =  definitions.rna_chi_pur
                if(mod.residues[k].res_mytype == "rU" or mod.residues[k].res_mytype == "rC"):
                    atoms = definitions.rna_chi_pyr

                cc_tmp = [mod.residues[k][at]  for at in atoms]
                coords_chi.append(cc_tmp)
                                    
            coords_bb = np.array(coords_bb)
            bb_torsion = bb_angles(coords_bb)

            coords_chi = np.array(coords_chi)
            chi_torsion = chi_angles(coords_chi)
            bb_torsion = np.degrees(np.column_stack((bb_torsion,chi_torsion)))

            # calculate angles
            #bb_torsion = cur_pdb.models[j].get_bb_torsions()
            if(args.hread):
                string = '# ' + files[i] + "-" + str(j) + "\n"
                string += "# " + "".join(cur_pdb.models[j].sequence) + "\n"
                for k in range(len(bb_torsion)):
                    string += "%10s" % (cur_pdb.models[j].sequence_id[k])
                    for l in range(len(bb_torsion[k])):
                        string += "%10.3f" % (bb_torsion[k][l])
                    string += "\n"
            else:
                string = files[i] + "." + str(j) + " "
                for k in range(len(bb_torsion)):
                    for l in range(len(bb_torsion[k])):
                        string += "%10.3f" % (bb_torsion[k][l])

                string += "\n"
            fh.write(string)
                


    fh.close()
    return 0
            
##################### ANNOTATE #######################
