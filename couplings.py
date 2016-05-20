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
import numpy as np
import mdtraj as md
import btools as bt
import definitions


def calc_j(angle,par):

    # standard Karplus
    cos = np.cos(angle+par[4])
    sin = np.sin(angle+par[4])
    val =  cos*cos*par[0] +  cos*par[1] + par[2]  + par[3]*cos*sin 
    #return angle
    return val

def stringify(angles,miss,seq,hread,raw):

    dd = angles.shape[1]/len(seq)
    result = []
    for t in xrange(angles.shape[0]):
        string = ""
        for ii in xrange(len(seq)):

            i1 = dd*ii
            aa =[angles[t,i1+jj] if ([jj,ii] not in miss) else float('nan') for jj in range(dd)]
            if(raw):
                stri = "".join(["%10.4f " %el for el in aa])
            else:
                aa1 = [calc_j(aa[jj],definitions.j3[jj][2]) for jj in range(dd)]
                stri = "".join(["%10.4f " %el for el in aa1])

                
            if(hread):
                stri = "%10s %s \n " % (seq[ii], stri)
            string += stri
        result.append(string)
        
    return result

def dihedral_idx(topology):

    miss = []
    indeces = []
    seq = []
    
    for rr in topology.residues:
        if(rr.name in definitions.rna):

            idxs = [[0,0,0,0] for x in range(len(definitions.j3))]
            seq.append("%s_%d_%d" % (rr.name,rr.resSeq,rr.chain.index))
            jj = len(seq)-1

            for j in range(len(definitions.j3)):
                r_atoms = definitions.j3[j][1]
                if(definitions.j3[j][0] == "C4Pe" or definitions.j3[j][0]=="H3P"):
                    try:
                        rr_p = topology.residue(rr.index+1)
                        idxs[j] = [rr.atom(r_atoms[0]).index,rr.atom(r_atoms[1]).index,\
                                rr.atom(r_atoms[2]).index,rr_p.atom(r_atoms[3]).index]
                    except:
                        miss.append([j,jj])
                else:
                    try:
                        idxs[j] = [rr.atom(r_atoms[0]).index,rr.atom(r_atoms[1]).index,\
                                rr.atom(r_atoms[2]).index,rr.atom(r_atoms[3]).index]
                    except:
                        miss.append([j,jj])
            indeces.extend(idxs)

    return indeces, seq, miss

def couplings(args):

    print "# Calculating Jcouplings angles..."
    fh = open(args.name,'a')

    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]
        
    for i in xrange(0,len(files)):
        print "#",files[i]

        if(args.pdbs!=None):
            cur_pdb = md.load_pdb(files[i])
            idxs,seq,miss = dihedral_idx(cur_pdb.topology)
            angles = md.compute_dihedrals(cur_pdb,idxs,periodic=False)
            string = stringify(angles,miss,seq,args.hread,args.raw)
            if(args.hread):
                fh.write("# file %s \n" % (files[i]))
                fh.write(string[0])
            else:
                fh.write("%s %s \n" % (files[i],string[0]))

        else:
            cur_pdb = md.load_frame(files[i],0,top=args.top)
            idxs,seq,miss = dihedral_idx(cur_pdb.topology)

            for chunk in md.iterload(files[i], chunk=100,top=args.top):
                
                angles = md.compute_dihedrals(chunk,idxs,periodic=False)
                string = stringify(angles,miss,seq,args.hread,args.raw)
            
                for t in range(len(string)):
                    if(args.hread):
                        fh.write("# file %s time=%10f \n" % (files[i], chunk[t].time))
                        fh.write(string[t])
                    else:
                        fh.write("%10f %s \n" % (chunk[t].time,string[t]))


        
            
    return 0
            

