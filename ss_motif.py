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

def ss_motif(args):

    # sanity checks
    files = args.files
    print "# Finding 3D Single Strand Motifs..."
    assert args.bulges < 3, "# FATAL: cannot do bulges > 2"

    ref_pdb = reader.Pdb(args.reference,res_mode=args.res_mode,at_mode="LCS")
    assert len(ref_pdb.models)==1, "# FATAL: Reference PDB file cannot have multiple models" 
    ref_len = len(ref_pdb.models[0].sequence)
    if(args.seq==None):
        query="N"*ref_len
    else:
        assert len(args.seq)==ref_len, "# FATAL: query structure and sequence length mismatch!"
        query = args.seq



    # OK...
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    ref_mat = ref_pdb.models[0].get_4dmat(args.cutoff)
    ref_mat_f = ref_mat.reshape(-1,4)

    treshold_sq = args.treshold*args.treshold 

    for i in xrange(0,len(files)):

        # if no PDBdump, read base atoms only (gives little speed-up)
        if(args.dump_pdb==True):
            cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode,at_mode="AA")
        else:
            cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode,at_mode="LCS")

        # read models
        counter = 0    
        for j in xrange(len(cur_pdb.models)):
            if(len(cur_pdb.models[j].sequence) < ref_len):
                continue
            
            all_idx = cur_pdb.models[j].get_idx(query,args.bulges)
            for idx in all_idx:
                gmat = cur_pdb.models[j].get_4dmat(args.cutoff,idx)
                red_mat = gmat.reshape(-1,4)
                diff = (red_mat-ref_mat_f)**2
                ermsd_sq = np.sum(np.sum(diff))/ref_len

                if(ermsd_sq < treshold_sq):
                    ermsd = np.sqrt(ermsd_sq)
                    seq = "_".join([cur_pdb.models[j].sequence_id[p] for p in idx ])
                    string = '%8.5f %s %i - %s \n' % (ermsd,files[i],j,seq)
                    fh.write(string)

                    if(args.dump_pdb == True):
                        counter += 1
                        seq = "".join([cur_pdb.models[j].sequence[p] for p in idx ])
                        ermsd_str = "%5.3f" % ermsd
                        new_pdb = args.name + "_" + files[i].split("/")[-1].split(".pdb")[0] + "_" + str(counter).zfill(4) + "_" + ermsd_str.strip() + ".pdb"
                        fh_pdb = open(new_pdb,'w')
                        fh_pdb.write(cur_pdb.models[j].string_pdb(idx))
                        fh_pdb.close()
        fh.close()
