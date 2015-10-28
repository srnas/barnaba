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
from scipy.spatial import distance
import tools

####################### MOTIF #########################

def ss_motif(args):

    # sanity checks
    files = args.files
    print "# Finding 3D Single Strand Motifs..."
    #assert args.bulges < 3, "# FATAL: cannot do bulges > 2"
    ref_pdb = reader.Pdb(args.reference,res_mode=args.res_mode)
    ref_len = len(ref_pdb.model.sequence)
    ref_mat = (ref_pdb.model.get_gmat(args.cutoff)).reshape(-1)
    if(args.seq==None):
        query="N"*ref_len
    else:
        assert len(args.seq)==ref_len, "# FATAL: query structure and sequence length mismatch!"
        query = args.seq

    pattern = tools.get_pattern(query)
    # OK...
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    for i in xrange(0,len(files)):

        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)
        cur_pdb.set_xtc(args.xtc)
        cur_len = len(cur_pdb.model.sequence)

        if(cur_len<ref_len): continue

        # return indeces matching query
        indeces = tools.get_idx(cur_pdb.model.sequence,query,args.bulges)
        idx = 0
        while(idx>=0):
            gmats = [(cur_pdb.model.get_gmat(args.cutoff,index)).reshape(-1) for index in indeces]
            dists = distance.cdist([ref_mat],gmats)/np.sqrt(ref_len)
        
            below_t = (dists<args.treshold).nonzero()[1]
            for ss in below_t:
                seq = "_".join([cur_pdb.model.sequence_id[p] for p in indeces[ss] ])
                string = '%8.5f %s %i - %s \n' % (dists[0,ss],files[i],idx,seq)
                fh.write(string)

            idx = cur_pdb.read()

    fh.close()
