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

#import reader as reader
import definitions
import mdtraj as md
import btools as bt
import numpy as np

##################### ANNOTATE #######################

def get_string(fname,time,pairs,annotation,seq_id,hread):

                    
    if(hread):
        string = "# %30s t=%10s \n" % (fname,time)
        ss = ["%10s %10s %2s \n" % (seq_id[pairs[j][0]],seq_id[pairs[j][1]],annotation[j]) for j in  xrange(len(pairs))]
        string = "".join(ss)
    else:
        #very expensive - use only for short sequences...
        merged = [str(el[0]) + "_" + str(el[1]) for el in pairs]
        string = "%f " % (time)
        for i1 in xrange(len(seq_id)):
            for i2 in xrange(i1+1,len(seq_id)):
                symbol = ".."
                lab = str(i1) + "_" + str(i2)
                if(lab in merged):
                    symbol = annotation[merged.index(lab)]
                string += "%2s " % symbol
        string += "\n"
    return string

    
def pymol_script(pdbname,seq,pairs,interactions):
    
    def color(itype):
        
        join = "".join(itype)

        if("WCc" in itype):
            return "firebrick"
        if("GUt" in itype):
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
        colors.append(color(ints))

    for i in range(len(colors)):
        fh.write("set cartoon_ring_color, %s , resi %d \n" %  ( colors[i], seq[i]+1))  
        fh.write("set cartoon_ladder_color, %s , resi %d \n" %  ( colors[i], seq[i]+1))  
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

def calculate(mat,angles,sequence):

    ll = len(sequence)
    # get upper triangular
    dist_sq = mat[:,:,0]**2 + mat[:,:,1]**2 + mat[:,:,2]**2
    # find pairs with cutoff < 1.6 (upper triangular only...)
    m_idx1 = np.where(np.triu(dist_sq,1)>0.0)
    m_idx2 = (m_idx1[1],m_idx1[0])
    rho1 = mat[m_idx1][:,0]**2 + mat[m_idx1][:,0]**2
    rho2 = mat[m_idx2][:,0]**2 + mat[m_idx2][:,0]**2
    z1 = mat[m_idx1][:,2]**2
    z2 = mat[m_idx2][:,2]**2
    
    openings = []
    closings = []
    pairs = []
    anno = []
    
    for i in range(len(m_idx1[0])):
        idx1=m_idx1[0][i]
        idx2=m_idx1[1][i]

        id1 = idx1,idx2
        id2 = idx2,idx1
        # stacking
        int_type = "XX"
        # criteria on |z| > 0.2

        if(z1[i]> 0.04 and z2[i] > 0.04):
            # criteria on angles (40 < x < 140)

            if(np.abs(angles[id1])>0.76):
        
                
                # criteria on rho < 2.5 AA 
                #print q1,q2
                if(rho1[i] < 0.0625 or rho2[i] < 0.0625): 
                    
                    if(mat[id1][2] > 0.02):
                        if(mat[id2][2] < -0.02):
                            int_type = ">>"
                        else:
                            int_type = "><"
                    else:
                        if(mat[id2][2] > 0.02):
                            int_type = "<<"
                        else:
                            int_type = "<>"
                            
        if(z1[i]< 0.04 and z2[i] < 0.04):
            int_type = ""
            edges = [np.arctan2(mat[id1][1],mat[id1][0]), np.arctan2(mat[id2][1],mat[id2][0])]
            for edge in edges:
                if(edge > definitions.theta1 and edge <= definitions.theta2):
                    int_type += "W"
                else:
                    if(edge > definitions.theta2 or edge <= definitions.theta3):
                        int_type += "H"
                    else:
                        int_type += "S"
            # watson-watson can be watson-crick
            if(int_type == "WW"):
                r1 = sequence[m_idx1[0][i]][0]
                r2 = sequence[m_idx1[1][i]][0]
                if(definitions.complementary[r1]==r2):
                    val = bt.wc_gaussian(10.*mat[id1])*bt.wc_gaussian(10.*mat[id2])
                    if(val>1.0e-08):
                        int_type = "WC"
                else:
                    if((r1=="G" and r2=="U") or (r1=="U" and r2=="G")):
                        int_type = "GU"

            # cis/trans
            if(angles[id1]<0):
                int_type += "c"
            else:
                int_type += "t"

            # for generating dot-bracket
            if(int_type == "WCc"):
                assert (idx1 not in openings)
                assert (idx1 not in openings)
                assert (idx2 not in closings)
                assert (idx2 not in closings)
                openings.append(idx1)
                closings.append(idx2)
                
        pairs.append([idx1,idx2])
        anno.append(int_type)
    dotbracket = dotbr(openings,closings,ll)
    return pairs,anno, dotbracket
    
def annotate(args):

    print "# Annotating RNA structures..."
    fh = open(args.name,'a')

    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]
        
    for i in xrange(0,len(files)):
        print "#",files[i]

        if(args.pdbs!=None):
            cur_pdb = md.load_pdb(files[i])
        else:
            top=args.top
            cur_pdb = md.load_frame(files[i],0,top=args.top)
            
        cur_idx = bt.get_lcs_idx(cur_pdb.topology)

        # get list of residues
        rna_residues = [cur_pdb.topology.atom(at).residue  for at in cur_idx[0]]
        
        # get sequence
        rna_seq = [definitions.residue_dict[rr.name]  for rr in rna_residues]
        rna_seq_id = ["%s_%d_%d" % (rr.name,rr.resSeq,rr.chain.index)  for rr in rna_residues]

        if(args.pdbs!=None):
            
            # coordinates for LCS calculations
            c1 = cur_pdb.xyz[0,cur_idx[0]]
            c2 = cur_pdb.xyz[0,cur_idx[1]]
            c3 = cur_pdb.xyz[0,cur_idx[2]]
            
            mat,angles = bt.get_mat_annotation(c1,c2,c3,cutoff=1.58)
            if(len(mat)==0): continue
            pairs, annotation, dotbracket = calculate(mat,angles,rna_seq)

            fh.write(get_string(files[i],i,pairs,annotation,rna_seq_id,args.hread))
            
            if(args.pymol == True):
                res_idxs_full = [rr.index  for rr in rna_residues]
                pymol_script(files[i],res_idxs_full, pairs,annotation)
                
        else:
            
            # analyze trajectory in chunks of 100
            for chunk in md.iterload(files[i], chunk=100,top=top):
                for j in range(len(chunk)):
                    c1 = chunk.xyz[j,cur_idx[0]]
                    c2 = chunk.xyz[j,cur_idx[1]]
                    c3 = chunk.xyz[j,cur_idx[2]]
                    mat,angles = bt.get_mat_annotation(c1,c2,c3,cutoff=1.58)
                    if(len(mat)==0): continue
                    pairs, annotation, dotbracket = calculate(mat,angles,rna_seq)
                    fh.write(get_string(files[i],chunk[j].time,pairs,annotation,rna_seq_id,args.hread))
            

    fh.close()
    return 0
            
##################### ANNOTATE #######################
