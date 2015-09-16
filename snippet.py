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
import definitions
import tools 
from scipy.spatial import distance


def snippet(args):
    
    files = args.files
    # check query sequence
    for item in args.seq:
        if(item not in reader.Names.known_abbrev):
            print "# FATAL Error. Symbol ", item, " not known. Use ACGU NYR"
            return 1

    query = args.seq.split("%")
    assert len(args.seq.split("%")) < 2 , "# Fatal error: max 1 strand"


    fh = open(args.name,'w')
    print "# SPLIT..."
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    ll = [len(el) for el in query]
    
    for i in xrange(0,len(files)):

        try:
            cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode,permissive=True)
        except:
            print "# SKIPPING", files[i]
            continue
        cur_len = len(cur_pdb.model.sequence)
        if(cur_len<sum(ll)): continue

        # single strand
        indeces = tools.get_idx(cur_pdb.model.sequence,query[0],bulges=0)

        # check chain consistency - remove 
        tools.chain_consistency(indeces,cur_pdb.model.sequence_id)

        if(len(query)==2):
            indeces2 = tools.get_idx(cur_pdb.model.sequence,query[1],bulges=0)
            tools.chain_consistency(indeces2,cur_pdb.model.sequence_id)
            # to be done....
            
        eof = True
        
        new_pdb_r = files[i].split("/")[-1].split(".pdb")[0] + "_"
        while(eof):
            
            for index in indeces:

                seq_out = "".join([cur_pdb.model.sequence[res] for res in index])
                f_res = cur_pdb.model.sequence_id[index[0]]
                l_res = cur_pdb.model.sequence_id[index[-1]]

                new_pdb = new_pdb_r + seq_out + "_" + f_res + "_" + l_res + ".pdb"
                fh_pdb = open(new_pdb,'w')
                fh_pdb.write(cur_pdb.model.string_pdb(index,noP=True))
                fh_pdb.close()

            eof = cur_pdb.read()


    fh.close()
    return 0


