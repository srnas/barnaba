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
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)
        
    # calculate interaction matrix of the reference structure
    ref_pdb = reader.Pdb(args.reference,res_mode=args.res_mode)
    ref_len = len(ref_pdb.model.sequence)
    ref_mat = ref_pdb.model.get_gmat(args.cutoff)

    if(args.xtc!=None):
        assert len(files)==1, "# Error: when providing XTC trajectories, specify a single reference PDB file with -f"
        
    #all the rest and calculate ERMSD on-the-fly
    for i in xrange(0,len(files)):
        
        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)
        cur_pdb.set_xtc(args.xtc)
        cur_len = len(cur_pdb.model.sequence)
        assert cur_len==ref_len, "# Fatal error: sequence lengths mismatch %d %d \n" %(cur_len,ref_len)

            
        idx = 0
        eof = True
        while(eof):
            cur_mat = cur_pdb.model.get_gmat(args.cutoff)
            diff = (ref_mat-cur_mat)**2
            val = np.sqrt(np.sum(diff)/ref_len)

            string = "%8i %10.6f  " % (idx,val)
            if(args.perres):
                string += ";"
                per_res = np.sqrt(np.sum(np.sum(diff,axis=2),axis=1)/ref_len)
                for k in xrange(len(per_res)):
                    string += " %10.6f " % (per_res[k])
            string += "%s \n" % (files[i])
            fh.write(string)
            
            idx += 1
            if(args.xtc==None):
                eof = cur_pdb.read()
            else:
                eof = cur_pdb.read_xtc()
                                
    fh.close()
    return 0


