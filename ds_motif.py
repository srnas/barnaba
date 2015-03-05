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


####################### MOTIF #########################

def ds_motif(args):

    # sanity checks
    files = args.files
    print "# Finding 3D Single Strand Motifs..."
    assert args.bulges < 3, "# FATAL: cannot do bulges > 2"

    ref_pdb = reader.Pdb(args.reference,base_only=True)
    assert len(ref_pdb.models)==1, "# FATAL: The query PDB contains more that one model"
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
    for k in args.__dict__:
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    ref_mat_tot = ref_pdb.models[0].get_4dmat(args.cutoff)
    ref_mat_tot_f = ref_mat_tot.reshape(-1,4)

    indeces1=np.arange(0,ref_len1)
    ref_mat1 = ref_pdb.models[0].get_4dmat(args.cutoff,indeces1)
    
    ref_mat_f1 = ref_mat1.reshape(-1,4)
    indeces2=np.arange(ref_len1,ref_len)
    ref_mat2 = ref_pdb.models[0].get_4dmat(args.cutoff,indeces2)
    
    ref_mat_f2 = ref_mat2.reshape(-1,4)

    # calculate center of mass distances
    # this will be used to prune the search!
    ref_com1 = np.sum(ref_pdb.models[0].get_com(indeces1),axis=0)/ref_len1
    ref_com2 = np.sum(ref_pdb.models[0].get_com(indeces2),axis=0)/ref_len2
    diff_com = (ref_com2-ref_com1)**2
    dd= np.sqrt(np.sum(diff_com))

    for i in xrange(0,len(files)):

        # if no PDBdump, read base atoms only (gives little speed-up)
        if(args.dump_pdb==True):
            cur_pdb = reader.Pdb(files[i],base_only=False)
        else:
            cur_pdb = reader.Pdb(files[i],base_only=True)

        # read models
        counter = 0    
        treshold_sq = args.treshold*args.treshold 

        for j in xrange(len(cur_pdb.models)):
            if(len(cur_pdb.models[j].sequence) < ref_len):
                continue

            # strand 1            
            tmp_idx1 = []
            com1 = []
            all_idx1 = cur_pdb.models[j].get_idx(query1,args.bulges)
            for idx in all_idx1:
                gmat = cur_pdb.models[j].get_4dmat(args.cutoff,idx,permissive=False)
                red_mat = gmat.reshape(-1,4)
                diff = (red_mat-ref_mat_f1)**2
                ermsd_sq = np.sum(np.sum(diff))/ref_len1
                if(ermsd_sq < treshold_sq):
                    tmp_idx1.append(idx)
                    com = np.sum(cur_pdb.models[j].get_com(idx),axis=0)/ref_len1
                    com1.append(com)

            # strand 2
            tmp_idx2 = []
            com2 = []
            all_idx2 = cur_pdb.models[j].get_idx(query2,args.bulges)
            for idx in all_idx2:
                gmat = cur_pdb.models[j].get_4dmat(args.cutoff,idx,permissive=False)
                red_mat = gmat.reshape(-1,4)
                diff = (red_mat-ref_mat_f2)**2
                ermsd_sq = np.sum(np.sum(diff))/ref_len2
                if(ermsd_sq < treshold_sq):
                    tmp_idx2.append(idx)
                    com = np.sum(cur_pdb.models[j].get_com(idx),axis=0)/ref_len2
                    com2.append(com)
                                        

             
            dmat = distance.cdist(com1,com2)
            c_idx = (dmat<2.5*dd).nonzero()
            for idx in range(len(c_idx[0])):
                idx_combo = tmp_idx1[c_idx[0][idx]]+tmp_idx2[c_idx[1][idx]]
                # skip overlapping
                if(len(np.union1d(tmp_idx1[c_idx[0][idx]],tmp_idx2[c_idx[1][idx]]))!=ref_len):
                    continue

                gmat = cur_pdb.models[j].get_4dmat(args.cutoff,idx_combo)
                red_mat = gmat.reshape(-1,4)
                diff = (red_mat-ref_mat_tot_f)**2
                ermsd_sq = np.sum(np.sum(diff))/ref_len
                if(ermsd_sq < treshold_sq):
                    ermsd = np.sqrt(ermsd_sq)
                    seq = "_".join([cur_pdb.models[j].sequence_id[p] for p in idx_combo ])
                    string = '%8.5f %s %i - %s \n' % (ermsd,files[i],j,seq)
                    fh.write(string)

                    #print  " ".join([cur_pdb.models[j].sequence_id[p] for p in idx_combo ]), np.sqrt(ermsd_sq)
                    if(args.dump_pdb == True):
                        counter += 1
                        seq = "".join([cur_pdb.models[j].sequence[p] for p in idx_combo ])
                        ermsd_str = "%5.3f" % ermsd
                        new_pdb = args.name + "_" + files[i].split("/")[-1].split(".pdb")[0] + "_" + str(counter).zfill(4) + "_" + ermsd_str.strip() + ".pdb"
                        fh_pdb = open(new_pdb,'w')
                        fh_pdb.write(cur_pdb.models[j].string_pdb(idx_combo))
                        fh_pdb.close()


    fh.close()


