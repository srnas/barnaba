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
import kde as kde
import tools as t

####################### SCORING ########################
        
def score(args,files):

    fh = open(args.name,'w')

    pb.write_args(args,fh)

    ref_atoms,ref_sequence = pb.get_coord(args.ff)
    ref_mat = []
    for ii,model in enumerate(ref_atoms):
        ref_lcs,ref_origo = t.coord2lcs(ref_atoms[ii])
        ref_mat_ii,ids = t.lcs2mat_score(ref_lcs,ref_origo,args.cutoff)
        ref_mat.extend(ref_mat_ii)

    kernel = kde.gaussian_kde(ref_mat)
    kernel.set_bandwidth(0.25)

    print "# KDE computed. Bandwidth=",kernel.factor
    print "# Calculating ESCORE..."

    for f in files:
        atoms,sequence = pb.get_coord(f)
        for ii,model in enumerate(atoms):
            #print model[0]

            lcs,origo = t.coord2lcs(model)
            # Cutoff is slightly augmented 
            #mat,ids = t.lcs2mat_score(lcs,origo,args.cutoff+0.2)
            mat,ids = t.lcs2mat_score(lcs,origo,args.cutoff+0.2)
            val = kernel(mat)
            #for kk in range(len(sequence[ii])):
            #    ss = 0.0
            #    for el in range(len(val)):
            #        if(kk in ids[el]):
            #            ss += val[el]
            #    print sequence[ii][kk],ss

            string = '%8.5f %s.%i \n' % (sum(val),f,ii)
            fh.write(string)

####################### SCORING ########################
