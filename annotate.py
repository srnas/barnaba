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


##################### ANNOTATE #######################

def pymol_script(pdbname,seq,interactions):
    
    def color(itype):

        join = "".join(itype)
        if("WC" in itype):
            return "firebrick"
        if("W" in join or "H" in join or "S" in join):
            return "yelloworange"
        if(">" in join or "<" in join):
            return "forest"
        return "blue"

    root = pdbname.split("/")[-1].split(".pdb")[0]

    fh = open(root + ".pml",'w')
    fh.write("load " + pdbname + ", tutto \n")
    fh.write("bg_color white\n")
    fh.write("hide everything \n")
    fh.write("show cartoon \n")
    fh.write("set orthoscopic, on \n")
    fh.write("set cartoon_color, gray60 \n")
    fh.write("set cartoon_nucleic_acid_color, gray60 \n")
    fh.write("cartoon oval \n")
    fh.write("set cartoon_ring_mode, 1 \n")
    fh.write("set cartoon_ring_finder, 2 \n")
    fh.write("set cartoon_ring_width, 0.3 \n")
    fh.write("set cartoon_ring_transparency, 0.2 \n")
    fh.write("set cartoon_oval_length, 0.4 \n")
    fh.write("set cartoon_oval_width, 0.2 \n")
    fh.write("set cartoon_nucleic_acid_mode \n")
    
    colors = []
    for res_idx in range(len(seq)):
        ints = []
        for ii in interactions:
            if(res_idx in ii):
                ints.append(ii[2])
        #print seq[res_idx],ints,
        colors.append(color(ints))
        #print colors[-1]
    for i in range(len(colors)):
        fh.write("set cartoon_ring_color, " + colors[i] + ", resi " + seq[i].split("_")[0] + " and chain \"" + seq[i].split("_")[2] + "\"\n")  
        fh.write("set cartoon_ladder_color, " + colors[i] + ", resi " + seq[i].split("_")[0] + " and chain \"" + seq[i].split("_")[2] + "\"\n")  


        #print "set cartoon_ring_color, "  + cc[p] + ", resi " + str(j+1+offset) + " and  struct_" + str(h)
    #print "set cartoon_ladder_color, " + cc[p] + ", resi " + str(j+1+offset) + " and  struct_" + str(h)
    fh.close()

def annotate(args):

    print "# Annotating RNA structures..."
    print "# Annotation is currently available in RNA-only mode..."

    files = args.files
    fh = open(args.name,'w')
    fh.write("# This is a baRNAba run.\n")
    for k in sorted(args.__dict__):
        s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
        fh.write(s)

        
    for i in xrange(0,len(files)):
        cur_pdb = reader.Pdb(files[i],res_mode="R",at_mode="LCS")

        for j in xrange(len(cur_pdb.models)):

            string = '# ' + files[i] + "-" + str(j) + "\n"
            
            # calculate interactions
            interactions,openings,closings = cur_pdb.models[j].get_annotation()

            # print sequence
            string += "# " + "".join(cur_pdb.models[j].sequence) + "\n"

            # check pseudoknots
            ll = len(cur_pdb.models[j].sequence)
            anno = ['.']*ll
            levels = [-1]*ll
            for idx1 in xrange(len(openings)):
                start1 = openings[idx1]
                end1 = closings[idx1]
                # up one level
                levels[start1] += 1
                levels[end1] += 1

                for idx2 in xrange(len(closings)):
                    end2 = closings[idx2]
                    if(levels[end2] == levels[start1]):
                        if(end2 > start1 and end2 < end1):
                            levels[start1] += 1
                            levels[end1] += 1

            for idx1 in xrange(len(openings)):
                start1 = openings[idx1]
                end1 = closings[idx1]

                anno[start1] = definitions.op[levels[start1]]
                anno[end1] = definitions.cl[levels[end1]]

            # print dotbracket
            string += '# ' + "".join(anno) + "\n"

            if(args.hread):
                # print interactions 
                for el in interactions:
                    string += "%10s %10s %2s \n" % (cur_pdb.models[j].sequence_id[el[0]],cur_pdb.models[j].sequence_id[el[1]],el[2])
                fh.write(string)
            else:
                merged = [str(el[0]) + "_" + str(el[1]) for el in interactions]
                for i1 in xrange(len(cur_pdb.models[j].sequence_id)):
                    for i2 in xrange(i1+1,len(cur_pdb.models[j].sequence_id)):
                        symbol = ".."
                        lab = str(i1) + "_" + str(i2)
                        if(lab in merged):
                            symbol = interactions[merged.index(lab)][2]
                        string += "%2s " % symbol
                fh.write(string + "\n")

            if(args.pymol == True):
                pymol_script(files[i],cur_pdb.models[j].sequence_id, interactions)

    fh.close()
    return 0
            
##################### ANNOTATE #######################
