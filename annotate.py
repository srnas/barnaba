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
import numpy as np

##################### ANNOTATE #######################

def pymol_script(pdbname,seq,pairs,interactions):
    
    def color(itype):

        join = "".join(itype)
        if("WC" in itype):
            return "firebrick"
        if("GU" in itype):
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
    for res_idx in xrange(len(seq)):
        ints = []
        for ii in xrange(len(pairs)):
            if(res_idx in pairs[ii]):
                ints.append(interactions[ii])
        #print seq[res_idx],ints,
        colors.append(color(ints))
        #print colors[-1]
    for i in range(len(colors)):
        fh.write("set cartoon_ring_color, " + colors[i] + ", resi " + seq[i].split("_")[0] + " and chain \"" + seq[i].split("_")[2] + "\"\n")  
        fh.write("set cartoon_ladder_color, " + colors[i] + ", resi " + seq[i].split("_")[0] + " and chain \"" + seq[i].split("_")[2] + "\"\n")  


        #print "set cartoon_ring_color, "  + cc[p] + ", resi " + str(j+1+offset) + " and  struct_" + str(h)
    #print "set cartoon_ladder_color, " + cc[p] + ", resi " + str(j+1+offset) + " and  struct_" + str(h)
    fh.close()


def dotbr(openings,closings,ll):

    # check pseudoknots
    dotbr = ['.']*ll
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
        
        dotbr[start1] = definitions.op[levels[start1]]
        dotbr[end1] = definitions.cl[levels[end1]]
    return "".join(dotbr)

def calculate(mat,idx,sequence):


    rho = np.sqrt(mat[:,0]**2 + mat[:,1]**2)
    dotp_scale = mat*np.array(definitions.scale)[np.newaxis,:]
    dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
    ll = len(sequence)

    # first, find pairs : d(ij)<cutoff, d(ji)< cutoff
    
    anno = []
    pairs = []
    for i1 in xrange(len(mat)):
        if(dotp_norm[i1] > 1.58):
            continue
        a1 = idx[i1][0]
        a2 = idx[i1][1]
        if(a2 > a1):
            pairs.append([a1,a2])
            anno.append([mat[i1]])
        else:
            try:
                ii = pairs.index([a2,a1])
                anno[ii].append(mat[i1])
            except:
                continue

    # now find annotation
    new_pairs = []
    new_anno = []

    # this is for dot-bracket
    openings = []
    closings = []
    for ii in xrange(len(anno)):
        coords =  anno[ii]
        if(len(coords) != 2):
            continue
        
        rho1 = coords[0][0]**2 + coords[0][1]**2
        rho2 = coords[1][0]**2 + coords[1][1]**2
        z1sq = coords[0][2]**2
        z2sq = coords[1][2]**2

        # stackings        
        if(z1sq>4 and z2sq>4):
            if(rho1<25 and rho2<25):
                if(coords[0][2] > 2.0):
                    if(coords[1][2] < -2.0):
                        int_type = ">>"
                    else:
                        int_type = "><"
                    new_pairs.append(pairs[ii])
                    new_anno.append(int_type)
                    continue
                else:
                    if(coords[1][2] > 2.0):
                        int_type = "<<"
                    else:
                        int_type = "<>"
                    new_pairs.append(pairs[ii])
                    new_anno.append(int_type)
                    continue
        # pairing
        else:

            angles = [np.arctan2(coords[0][1],coords[0][0]),np.arctan2(coords[1][1],coords[1][0])]
            int_type = ""
            for angle in angles:
                if(angle > definitions.theta1 and angle <= definitions.theta2):
                    int_type += "W"
                else:
                    if(angle > definitions.theta2 or angle <= definitions.theta3):
                        int_type += "H"
                    else:
                        int_type += "S"
                        
            # watson-watson can be watson-crick
            if(int_type == "WW"):
                r1 = sequence[pairs[ii][0]]
                r2 = sequence[pairs[ii][1]]
                if(tools.is_complementary(r1,r2)):
                    val = tools.wc_gaussian(coords[0])*tools.wc_gaussian(coords[1])
                    if(val>1.0e-08):
                        int_type = "WC"
                else:
                    if((r1=="G" and r2=="U") or (r1=="U" and r2=="G")):
                        if(np.abs(coords[1][2]) <2.0 and np.abs(coords[0][2]) <2.0):
                            int_type = "GU"
                if(int_type == "WC"):
                    assert (pairs[ii][0] not in openings)
                    assert (pairs[ii][1] not in openings)
                    assert (pairs[ii][0] not in closings)
                    assert (pairs[ii][1] not in closings)
                    openings.append(pairs[ii][0])
                    closings.append(pairs[ii][1])
                  
            new_pairs.append(pairs[ii])
            new_anno.append(int_type)
    dotbracket = dotbr(openings,closings,ll)
    return new_pairs, new_anno, dotbracket
    
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
        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)

        if(args.xtc!=None):
            cur_pdb.set_xtc(args.xtc)
        idx = 0
        eof = True
        seq_id = cur_pdb.model.sequence_id
        string = "#" + files[i]  + "\n"
        string += "# " + "".join(cur_pdb.model.sequence) + "\n"

        while(eof):

            cutoff = 1.58
            mat,m_idx = cur_pdb.model.get_3dmat(cutoff,range(len(cur_pdb.model.sequence)))
            pairs,annotation,dotbracket = calculate(mat,m_idx,cur_pdb.model.sequence)

            string += "# " + dotbracket + "\n"
            if(args.hread):
                
                # print interactions
                
                for j in xrange(len(pairs)):
                    string += "%10s %10s %2s \n" % (seq_id[pairs[j][0]],seq_id[pairs[j][1]],annotation[j])
                fh.write(string)
            else:
                # very expensive - use only for short sequences...
                merged = [str(el[0]) + "_" + str(el[1]) for el in pairs]
                for i1 in xrange(len(seq_id)):
                    for i2 in xrange(i1+1,len(seq_id)):
                        symbol = ".."
                        lab = str(i1) + "_" + str(i2)
                        if(lab in merged):
                            symbol = annotation[merged.index(lab)]
                        string += "%2s " % symbol
                fh.write(string + "\n")

            if(args.pymol == True):
                pymol_script(files[i],seq_id, pairs,annotation)

            idx += 1
            if(args.xtc==None):
                eof = cur_pdb.read()
            else:
                eof = cur_pdb.read_xtc()

    fh.close()
    return 0
            
##################### ANNOTATE #######################
