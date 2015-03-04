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


def ermsd(args):
    
    files = args.files
    print "# Calculating ERMSD..."

    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in args.__dict__:
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)


    # calculate interaction matrix of the reference structure
    ref_pdb = reader.Pdb(args.reference,base_only=True)
    assert len(ref_pdb.models)==1, "# FATAL: Reference PDB file cannot have multiple models" 
    ref_len = len(ref_pdb.models[0].sequence)
    ref_mat = ref_pdb.models[0].get_4dmat(args.cutoff)
    ref_mat_f = ref_mat.reshape(-1,4)

    #all the rest and calculate ERMSD on-the-fly
    for i in xrange(0,len(files)):
        cur_pdb = reader.Pdb(files[i],base_only=True)
        for j in xrange(len(cur_pdb.models)):
            assert ref_len == len(cur_pdb.models[j].sequence) , "# FATAL: sequence lenghts mismatch"

            # get 4d matrix
            cur_mat = cur_pdb.models[j].get_4dmat(args.cutoff)
            cur_mat_f = cur_mat.reshape(-1,4)
            
            # calculate difference
            diff = (cur_mat_f-ref_mat_f)**2
            ermsd = np.sqrt(np.sum(np.sum(diff))/ref_len)

            if(args.ermsf==False):
                string = '%10.6f %s - %i \n' % (ermsd,files[i],j)
            else:
                string = '%10.6f - ' % (ermsd)
                for k in xrange(ref_len):
                    diff1 = (cur_mat[k,:]-ref_mat[k,:])**2
                    ermsf = np.sqrt(np.sum(np.sum(diff1))/ref_len) 
                    string += " %10.6f " % (ermsf)
                string += '- %s - %i \n' % (files[i],j)
            fh.write(string)
        print "# Calculated", len(cur_pdb.models),"ERMSD",files[i]

    fh.close()
    return 0


