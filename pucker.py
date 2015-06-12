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
import sys

def pucker_angles(coords):

    # duplicate some coordinates to allow looping
    coords_plus = np.hstack((coords,coords[:,0:3]))
    diffs = coords_plus[:,:-1]-coords_plus[:,1:]
    norm_sq = np.sum(diffs**2,axis=2)
    # if one atom is missing, all are set to nan
    idx_no = np.where(norm_sq>definitions.maxbond_sq)[0]
    norm_sq[idx_no,:] = 0.0

    # calculate angles. The mathc could be in principle moved to a
    # generic calc_dihedral function... however with this approach 
    # norm is calculated only once

    angles = []
    for i in xrange(coords_plus.shape[1]-3):
        b0 = -diffs[:,i]
        b1 = diffs[:,i+1]
        b2 = diffs[:,i+2]
        ss = np.sum(b0*b1,axis=1)

        v0= b0-((np.sum(b0*b1,axis=1)*b1.T)/norm_sq[:,i+1]).T
        v2= b2-((np.sum(b2*b1,axis=1)*b1.T)/norm_sq[:,i+1]).T
        x = np.sum(v0*v2,axis=1)
        m = np.cross(v0,b1).T/np.sqrt(norm_sq[:,i+1])
        y = np.sum(m.T*v2,axis=1)
        angles.append(np.arctan2( y, x ))
    angles = np.array(angles)

    # now calculate phase and amplitude of pucker
    x1 = angles[4]+angles[1]-angles[3]-angles[0]
    #Pconst = 2.0*(np.sin(np.pi/5.) + np.sin(np.pi/2.5))
    x2 = 3.0776835*angles[2]

    # phase
    p0 = np.arctan2(x1,x2) 
    p0[np.where( p0 < 0.0 )] += 2.0*np.pi

    # amplitude
    tm = angles[2]/np.cos(p0)  

    #print angles.shape,p0.shape,tm.shape
    angles = np.vstack((angles,tm))
    angles = np.vstack((angles,p0))

    return np.degrees(angles.T)


    
def pucker(args):

    print "# Sugar pucker ..."
    files = args.files
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    for i in xrange(0,len(files)):

        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode,at_mode="PUCKER")

        for j in xrange(len(cur_pdb.models)):
           
            mod = cur_pdb.models[j]

            coords = []
            for k in xrange(len(mod.residues)):               
                cc = [mod.residues[k][at] for at in definitions.rna_pucker]
                coords.append(cc)
            coords = np.array(coords)
            pucker = pucker_angles(coords)


            if(args.hread):
                string = '# ' + files[i] + "-" + str(j) + "\n"
                string += "# " + "".join(cur_pdb.models[j].sequence) + "\n"
                for k in xrange(len(pucker)):
                    string += "%10s" % (cur_pdb.models[j].sequence_id[k])
                    for l in xrange(len(pucker[k])):
                        string += "%10.3f " % pucker[k][l]
                    string += "\n"

            else:

                string = files[i] + "." + str(j) + " "
                for k in xrange(len(pucker)):
                    for l in xrange(len(pucker[k])):
                        string += "%10.3f " % pucker[k][l]
                string += "\n"

            fh.write(string)
                
            # calculate interactions
            #ss_torsion = cur_pdb.models[j].get_sugar_torsion()



    fh.close()
    return 0
            
##################### ANNOTATE #######################
