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

import pdbreader as pb
import tools as t
import numpy as N

op = ['(','[','{','<']
cl = [')',']','}','>']          

##################### ANNOTATE #######################

def annotate(args,files):

    print "# Annotating RNA structures..."

    fh = open(args.name,'w')
    pb.write_args(args,fh)

    for f in files:
        atoms,sequence = pb.get_coord(f)
        count = 0
        for model,seq in zip(atoms,sequence):
            count += 1
            string = '# ' + f + "-" + str(count) + "\n"
            
            lcs,origo = t.coord2lcs(model)
            mat = t.lcs2mat_3d(lcs,origo,args.cutoff)
            int_mat = t.analyze_mat(mat,seq)

            anno_string = ''
            seq_string = '# SEQ '

            ll = len(seq)
            anno = ['.']*ll
            levels = [-1]*ll
            openings = []
            closings = []

            for j in range(int_mat.shape[0]):
                seq_string += seq[j] + " "
                for k in range(j+1,int_mat.shape[0]):
                    if(args.compact==True):
                        anno_string += "%4s " % (t.interactions[int(int_mat[j,k])])
                    else:
                        if(int_mat[j,k] != 0):
                            tt = t.interactions[int(int_mat[j,k])]
                            anno_string += "%10s %10s %4s \n" % (seq[j],seq[k],t.interactions[int(int_mat[j,k])])
                            r1 = seq[j].split("_")[1]
                            r2 = seq[k].split("_")[1]
                            if(tt == 'WC'):
                                if(j in openings or j in closings):
                                    print "# Warning - Residue", seq[j],"has double WC contact. This should not happen!"
                                if(k in openings or k in closings):
                                    print "# Warning - Residue", seq[k],"has double WC contact. This should not happen!"
                                openings.append(j)
                                closings.append(k)
                                
            # pseudoknots check
            for i in range(ll):

                if(i in openings):
                    idx = openings.index(i)
                    levels[i] += 1
                    anno[i] = op[levels[i]]
                    levels[closings[idx]] += 1
                    anno[closings[idx]] = cl[levels[i]]
                    for j in range(i+1,closings[idx]):
                        if(j in closings):
                            jdx = closings.index(j)
                            if(openings[jdx] < i and anno[j] == cl[levels[i]]):
                                levels[i] += 1
                                anno[i] = op[levels[i]]
                                levels[closings[idx]] += 1
                                anno[closings[idx]] = cl[levels[i]]
                                
            dotb_string = '# '
            for el in anno:
                dotb_string += el


            fh.write(string)
            fh.write(seq_string + "\n")
            fh.write(dotb_string + "\n")
            fh.write(anno_string)

    fh.close()
    return 0
            
##################### ANNOTATE #######################
