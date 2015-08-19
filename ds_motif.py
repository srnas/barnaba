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
import tools
from scipy.spatial import distance


####################### MOTIF #########################

def ds_motif(args):

    # sanity checks
    files = args.files
    print "# Finding 3D Double Stranded Motifs..."

    ref_pdb = reader.Pdb(args.reference,res_mode=args.res_mode)

    ref_len1 = args.l1
    ref_len2 = args.l2 
    ref_len = ref_len1+ref_len2
    if(args.seq==None):
        query1="N"*ref_len1
        query2="N"*ref_len2
    else:
        query1=(args.seq).split("%")[0]
        query2=(args.seq).split("%")[1]
        assert len(query1)==ref_len1, "# FATAL: query structure and sequence length mismatch!"
        assert len(query2)==ref_len2, "# FATAL: query structure and sequence length mismatch!"



    # OK...
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    ref_mat_tot = ref_pdb.model.get_gmat(args.cutoff).reshape(-1)

    indeces1=np.arange(0,ref_len1)
    ref_mat1 = ref_pdb.model.get_gmat(args.cutoff,indeces1).reshape(-1)
    
    indeces2=np.arange(ref_len1,ref_len)
    ref_mat2 = ref_pdb.model.get_gmat(args.cutoff,indeces2).reshape(-1)
    

    # calculate center of mass distances
    # this will be used to prune the search!
    ref_com1 = ref_pdb.model.get_lcs_com(indeces1)
    ref_com2 = ref_pdb.model.get_lcs_com(indeces2)
    diff_com = (ref_com2-ref_com1)**2
    dd= np.sqrt(np.sum(diff_com))

    for i in xrange(0,len(files)):

        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)
        cur_pdb.set_xtc(args.xtc)

        cur_len = len(cur_pdb.model.sequence)

        if(cur_len<ref_len): continue

        all_idx1 = tools.get_idx(cur_pdb.model.sequence,query1,args.bulges)
        all_idx2 = tools.get_idx(cur_pdb.model.sequence,query2,args.bulges)


        idx = 0
        eof = True
        while(eof):

            # get indeces and coms of first half
            gmats1 =  [(cur_pdb.model.get_gmat(args.cutoff,index)).reshape(-1) for index in all_idx1]
            dists1 = distance.cdist([ref_mat1],gmats1)/np.sqrt(ref_len1)
            below_t1 = (dists1<args.treshold).nonzero()[1]
            idx1 = [all_idx1[j] for j in below_t1]
            com1 = [cur_pdb.model.get_lcs_com(all_idx1[j]) for j in below_t1]
            
            # get indeces of second half
            gmats2 =  [(cur_pdb.model.get_gmat(args.cutoff,index)).reshape(-1) for index in all_idx2]
            dists2 = distance.cdist([ref_mat2],gmats2)/np.sqrt(ref_len2)
            below_t2 = (dists2<args.treshold).nonzero()[1]
            idx2 = [all_idx2[j] for j in below_t2]
            com2 = [cur_pdb.model.get_lcs_com(all_idx2[j]) for j in below_t2]
            
            # calculate all distances between center of mass
            dmat = distance.cdist(com1,com2)
            c_idx = (dmat<1.5*dd).nonzero()

            # get combo indeces
            dmine = [dmat[c_idx[0][ii],c_idx[1][ii]] for ii in range(len(c_idx[0]))]
            idx_combo = [idx1[c_idx[0][ii]] + idx2[c_idx[1][ii]] for ii in range(len(c_idx[0]))]
            gmatsf = [(cur_pdb.model.get_gmat(args.cutoff,index)).reshape(-1) for index in idx_combo]
            distsf = distance.cdist([ref_mat_tot],gmatsf)/np.sqrt(ref_len)
            below_tf = (distsf<args.treshold).nonzero()[1]
            for ss in below_tf:
                print dd,dmine[ss],1.5*dd
                seq = "_".join([cur_pdb.model.sequence_id[p] for p in idx_combo[ss] ])
                string = '%8.5f %s %i - %s \n' % (distsf[0,ss],files[i],idx,seq)
                fh.write(string)

            idx += 1
            
            if(args.xtc==None):
                eof = cur_pdb.read()
            else:
                eof = cur_pdb.read_xtc()



    fh.close()


