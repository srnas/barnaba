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
import kde as kde
#import tools as tools

####################### SCORING ########################
        
def score(args):

    files = args.files
    print "# Calculating ESCORE..."
    print "# This is a baRNAba run."

    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)


    # calculate interaction matrix of the reference structure
    ref_pdb = reader.Pdb(args.ff,res_mode=args.res_mode,at_mode="LCS")
    ref_mat = []
    for j in xrange(len(ref_pdb.models)):
        ref_mat_ii = ref_pdb.models[j].get_mat_score(args.cutoff)
        ref_mat.extend(ref_mat_ii)
    kernel = kde.gaussian_kde(ref_mat)
    kernel.set_bandwidth(0.25)

    print "# KDE computed. Bandwidth=",kernel.factor
    print "# Calculating ESCORE..."

    for i in xrange(0,len(files)):
        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode,at_mode="LCS")
        for j in xrange(len(cur_pdb.models)):
            mat = cur_pdb.models[j].get_mat_score(args.cutoff+0.2)
            val = kernel(mat)
            string = '%8.5f ' % (sum(val))
            string += '%s.%i \n' % (files[i],j)
            fh.write(string)

####################### SCORING ########################
