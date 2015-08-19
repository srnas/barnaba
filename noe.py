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

import numpy as np
import reader as reader

####################### MOTIF #########################

def noe(args):

    # sanity checks
    files = args.files
    print "# Calculating NOES"
    
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    for i in xrange(0,len(files)):

        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)
        cur_pdb.set_xtc(args.xtc)
        cur_pdb.model.set_h_idx()
        
        idx = 0
        eof = True
        data = []
        while(eof):
            data.append(cur_pdb.model.calc_pairwise_h())
            idx += 1
            if(args.xtc==None):
                eof = cur_pdb.read()
            else:
                eof = cur_pdb.read_xtc()

        data = np.array(data)
        triu = np.triu_indices(len(cur_pdb.model.h_labels), 1)
        string = ""
        for j in range(len(data[0])):
            avg = np.power(np.average(data[:,j]),-1./6.)
            if(avg<args.cutoff):
                string += "%12s %12s %8.4f \n" % (cur_pdb.model.h_labels[triu[0][j]],cur_pdb.model.h_labels[triu[1][j]], avg)
        fh.write(string)

    fh.close()
