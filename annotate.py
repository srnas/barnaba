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
            string = '# ' + f + "-" + str(count) + " "
            lcs,origo = t.coord2lcs(model)
            mat = t.lcs2mat_3d(lcs,origo,args.cutoff)
            int_mat = t.analyze_mat(mat,seq)
            if(args.compact==True):
                for j in range(int_mat.shape[0]):
                    for k in range(j+1,int_mat.shape[0]):
                        string += "%4s " % (t.interactions[int(int_mat[j,k])])
                fh.write(string + "\n")
            else:
                string += "\n"
                fh.write(string)
                for j in range(int_mat.shape[0]):
                    for k in range(j+1,int_mat.shape[0]):
                        if(int_mat[j,k] != 0):
                            string  = "%10s %10s %4s \n" % (seq[j],seq[k],t.interactions[int(int_mat[j,k])])
                            fh.write(string)
    fh.close()
    return 0
            
##################### ANNOTATE #######################
