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
#import tools as tools
from scipy.spatial import distance


def snippet(args):
    
    files = args.files
    # check query sequence
    for item in args.seq:
        if(item not in reader.Names.known_abbrev):
            print "# FATAL Error. Symbol ", item, " not known. Use ACGU NYR %"
            return 1

    query = args.seq.split("%")
    if(len(query)>2):
        print "# FATAL Error. Query sequence cannot be composed by more than 2 subsequences"
        return 1
        

    fh = open(args.name,'w')
    print "# SPLIT..."
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)


    for i in xrange(0,len(files)):
        pdb = reader.Pdb(files[i],res_mode=args.res_mode,at_mode="AA")
        for j in xrange(len(pdb.models)):
            # len 1 - easy!
            if(len(query)==1):
                if(len(pdb.models[j].sequence) < len(query[0])):
                    continue
                all_idx = pdb.models[j].get_idx(query[0])

            # len 2 - difficult!
            if(len(query)==2):
                if(len(pdb.models[j].sequence) < len(query[0]) + len(query[1])):
                    continue

                all_idx_1 = pdb.models[j].get_idx(query[0])
                all_idx_2 = pdb.models[j].get_idx(query[1]) 
                all_idx = []
                for el1 in all_idx_1:
                    for el2 in all_idx_2:

                        # check sequence overlap
                        if(el1[-1] < el2[0] or el2[-1] < el1[0]):
                            # check 3d distance
                            com1 = pdb.models[j].get_com(el1)
                            com2 = pdb.models[j].get_com(el2)
                            distances = distance.cdist(com1,com2)

                            # print if com distance is less than 15 A
                            if(distances.min()<15.0):  
                                all_idx.append(el1+el2)

                                #if(args.helix):
                                #    # check annottion
                                #    #print pdb.models[j].sequence_id[el1[0]], pdb.models[j].sequence_id[el2[0]]
                                #    pp,op,cl = pdb.models[j].get_annotation(el1+el2)
                                #    pair1= [0,len(el1)+len(el2)-1,"WC"]
                                #    #pair2= [len(el1),len(el1)+len(el2),"WC"]
                            # 
                            #        if(pair1 in pp):
                                #        all_idx.append(el1+el2)
                                #        #print pp, pair1," ".join([pdb.models[j].sequence_id[res] for res in el1+el2])
                                #                                else:
                                 

            for idx in all_idx:

                # print to output
                out = files[i] + "_model_" + str(j+1) + "_seq_" 
                out += "/".join([pdb.models[j].sequence_id[res] for res in idx])
                seq_out = "".join([pdb.models[j].sequence[res] for res in idx])



                if(args.dump_pdb == True):
                    new_pdb = files[i].split("/")[-1].split(".pdb")[0] + "_" + seq_out + "_" \
                        + pdb.models[j].sequence_id[idx[0]][:-2] + "_"\
                        + pdb.models[j].sequence_id[idx[-1]][:-2] + ".pdb"
                    fh_pdb = open(new_pdb,'w')
                    fh_pdb.write(pdb.models[j].string_pdb(idx))
                    fh_pdb.close()

                gmat = pdb.models[j].get_4dmat(args.cutoff,idx)
                fh.write(out + " ")            
                gout = ''
                for i1 in xrange(gmat.shape[0]):
                    for i2 in xrange(gmat.shape[1]):
                        if(i1!=i2):
                            for i3 in xrange(gmat.shape[2]):
                                gout += '%10.6f ' % gmat[i1,i2,i3]
                fh.write(gout + "\n")
                    
    fh.close()
    return 0


