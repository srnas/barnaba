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

    
#
#def bb_angles(coords):
#
#    diffs = coords_plus[:,:-1]-coords_plus[:,1:]
#    norm_sq = np.sum(diffs**2,axis=2)
#    #idx_up = np.where(norm_sq>definitions.maxbond_pucker_sq)
#    #idx_low = np.where(norm_sq<definitions.minbond_pucker_sq)
#    # calculate angles. The mathc could be in principle moved to a
#    # generic calc_dihedral function... however with this approach 
#    # norm is calculated only once
#
#    angles = []
#    for i in xrange(coords_plus.shape[1]-3):
#        b0 = -diffs[:,i]
#        b1 = diffs[:,i+1]
#        b2 = diffs[:,i+2]
#        ss = np.sum(b0*b1,axis=1)
#
#        v0= b0-((np.sum(b0*b1,axis=1)*b1.T)/norm_sq[:,i+1]).T
#        v2= b2-((np.sum(b2*b1,axis=1)*b1.T)/norm_sq[:,i+1]).T
#        x = np.sum(v0*v2,axis=1)
#        m = np.cross(v0,b1).T/np.sqrt(norm_sq[:,i+1])
#        y = np.sum(m.T*v2,axis=1)
#        angles.append(np.arctan2( y, x ))
#    angles = np.array(angles)
#
#    return np.degrees(angles.T),idx_up,idx_low
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
            coords_chi = np.array(coords_chi)
            
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
