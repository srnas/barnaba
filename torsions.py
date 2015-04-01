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

def torsions(args):

    print "# Annotating RNA structures..."
    files = args.files
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in args.__dict__:
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

    for i in xrange(0,len(files)):
        cur_pdb = reader.Pdb(files[i],base_only=False)
        for j in xrange(len(cur_pdb.models)):

            string = '# ' + files[i] + "-" + str(j) + "\n"
            
            # print sequence
            string += "# " + "".join(cur_pdb.models[j].sequence) + "\n"
            
            # calculate interactions
            bb_torsion = cur_pdb.models[j].get_bb_torsions()
            if(args.hread):
                for k in range(len(bb_torsion)):
                    string += "%10s" % (cur_pdb.models[j].sequence_id[k])
                    for l in range(len(bb_torsion[k])):
                        try:
                            string += "%10.3f" % bb_torsion[k][l]
                        except:
                            string += "%10s" % bb_torsion[k][l]
                    string += "\n"
            else:
                for k in range(len(bb_torsion)):
                    for l in range(len(bb_torsion[k])):
                        try:
                            string += "%10.3f" % bb_torsion[k][l]
                        except:
                            string += "%10s" % bb_torsion[k][l]
            string += "\n"
            fh.write(string)
                
            # calculate interactions
            #ss_torsion = cur_pdb.models[j].get_sugar_torsion()



    fh.close()
    return 0
            
##################### ANNOTATE #######################
