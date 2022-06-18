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

from __future__ import absolute_import, division, print_function
import sys
import glob
import os
import argparse
import barnaba as bb
import itertools as its

if (sys.version_info > (3, 0)):
   _HAS_PYTHON3=True
else:
   _HAS_PYTHON3=False

def parse(arguments):

    parser = argparse.ArgumentParser(description='This is baRNAba')
    subparsers = parser.add_subparsers(title="Subcommands",dest='subparser_name')
    subparsers.required = True
    
    # ERMSD PARSER
    parser_01 = subparsers.add_parser('ERMSD', help='calculate eRMSD')
    parser_01.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_01.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",required=False,default=None)
    parser_01.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_01.add_argument("--top", dest="top",help="Topology file",required=False)
    
    parser_01.add_argument("--ref", dest="reference",help="Reference PDB file",required=True)
    parser_01.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)

    
    # RMSD parser
    parser_01a = subparsers.add_parser('RMSD', help='calculate RMSD')
    parser_01a.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_01a.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",required=False,default=None)
    parser_01a.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_01a.add_argument("--top", dest="top",help="Topology file",required=False)    
    parser_01a.add_argument("--ref", dest="reference",help="Reference PDB file",required=True)
    parser_01a.add_argument("--heavy-atom", dest="heavy_atom",help="Attempt to use all heavy atoms for superposition",action='store_true',default=False)
    parser_01a.add_argument("--dump", dest="dump",help="Write aligned PDB/TRJ",action='store_true',default=False)

    
    # ESCORE PARSER 
    parser_02 = subparsers.add_parser('ESCORE', help='Calculate eSCORE')
    parser_02.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_02.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",required=False,default=None)
    parser_02.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_02.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_02.add_argument("--ff", dest="reference",help="Force-field PDB file",required=True)
    parser_02.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=1.58)",default=1.58,type=float)

    # SEARCH SINGLE STRANDED MOTIFS - OK
    parser_03 = subparsers.add_parser('SS_MOTIF', help='Search single stranded (hairpin loop) RNA motif ')
    parser_03.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_03.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_03.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_03.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_03.add_argument("--query", dest="query",help="Query PDB file",required=True)
    parser_03.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_03.add_argument("--threshold", dest="threshold",help="ERMSD threshold",default=0.7,type=float)    
    parser_03.add_argument("--bulges", dest="bulges",help="Number of allowed bulged nucleotides",default=0,type=int)   
    parser_03.add_argument("--sequence", dest="seq",help="Sequence Accepts ACGU/NRY/ format. Default = any",required=False,default=None)
    parser_03.add_argument("--dump", dest="dump",help="Write pdb files",action='store_true',default=False)

    # SEARCH DOUBLE STRANDED MOTIFS
    parser_04 = subparsers.add_parser('DS_MOTIF', help='Search double stranded RNA motifs ')
    parser_04.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_04.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_04.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_04.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_04.add_argument("--query", dest="query",help="Reference PDB file",required=True)
    parser_04.add_argument("--l1", dest="l1",help="Length of first strand",required=True,type=int)    
    parser_04.add_argument("--l2", dest="l2",help="Lenght of second strand",required=True,type=int)    

    parser_04.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_04.add_argument("--threshold", dest="threshold",help="ERMSD threshold",default=0.7,type=float)    
    parser_04.add_argument("--bulges", dest="bulges",help="Number of allowed bulged nucleotides",default=0,type=int)   
    parser_04.add_argument("--sequence", dest="seq",help="Sequence Accepts ACGU/NRY/ format. Default = any",required=False,default=None)
    parser_04.add_argument("--dump", dest="dump",help="Write pdb files",action='store_true',default=False)
    
    # ANNOTATE
    parser_05 = subparsers.add_parser('ANNOTATE', help='Annotate RNA structure')
    parser_05.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_05.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_05.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_05.add_argument("--top", dest="top",help="Topology file",required=False)
    parser_05.add_argument("--dotbracket", dest="dotbr",help="write dot-bracket annotation",action='store_true',default=False)


    # SEC_STRUCTURE
    parser_11 = subparsers.add_parser('SEC_STRUCTURE', help='Draw secondary structure from annotation')
    parser_11.add_argument("-o", dest="name", help="output_name",default=None,required=False)
    parser_11.add_argument("--ann", dest="f_anns", help="Annotation file(s) ([pairing and/or stacking] or dotbracket). Use (pairing and stacking) or dotbracket to minimize and show. Use all 3 files to minimize on dotbracket but show all annotations.", nargs="+",default=None,required=True)
    parser_11.add_argument("--draw_interm", dest="draw_interm", help="Draw intermediate structures (int)", default=0)
    parser_11.add_argument("--template", dest="template", help="SVG template structure (VARNA, RNAstructure etc.) ", default=None, required=False)
#    parser_11.add_argument("-T", dest="T_init", help="Initial temperature for annealing", default=0)
    parser_11.add_argument("--nsteps", dest="nsteps", help="Number of steps for geometry optimization", default=1000)
#    parser_11.add_argument("--no_tertiary_contacts", dest="tertiary_contacts", help="Do not use tertiary contacts in geometry optimization", action='store_false', default=True)
    parser_11.add_argument("--output_ids", dest="output_ids", help="Print residue IDs instead of letters",action='store_true',default=False)
#    parser_11.add_argument("--chi", dest="chi_conf", help="Chi conformation (int) translated from binary system (chi in syn =1, not syn =0) for all residues",default=0)
  #  parser_11.add_argument("-w", dest="f_weights", help="Weights file (one weight per frame)",default=None, required=False)
    parser_11.add_argument("--threshold", dest="thresh", help="Threshold for interaction frequency (default 0.1)",default=0.1, required=False)
    parser_11.add_argument("--not_order", dest="not_order", help="Don't try to make double strands rectangular", action='store_true',default=False)

    
    # DUMP
    parser_06 = subparsers.add_parser('DUMP', help='DUMP structural parameters')
    parser_06.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_06.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_06.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_06.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_06.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_06.add_argument("--dumpG", dest="dumpG",help="Write G vectors on .gvec file",action='store_true',default=False)
    parser_06.add_argument("--dumpR", dest="dumpR",help="Write R vectors on .rvec file",action='store_true',default=False)
    

    # CALCULATE TORSION ANGLES
    parser_08 = subparsers.add_parser('TORSION', help='Calculate dihedral angles')
    parser_08.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_08.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_08.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_08.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_08.add_argument("--backbone", dest="backbone",help="Calculate backbone (a,b,g,d,e,z) and chi torsion angle",action='store_true',default=False)
    parser_08.add_argument("--sugar", dest="sugar",help="Calculate sugar torsion angles (v0...v5)",action='store_true',default=False)
    parser_08.add_argument("--pucker", dest="pucker",help="Calculate sugar pucker pseudorotation parameters (amplitude, phase)",action='store_true',default=False)
    parser_08.add_argument("--altona", dest="altona",help="Calculate sugar pucker angles using Altona-Sundaralingam treatment",action='store_true',default=False)
    parser_08.add_argument("--res", dest="res",help="Calculate torsion angle for specific residues",required=False,nargs="+",default=None)
    

    # CALCULATE TORSION ANGLES - J COUPLINGS
    parser_09 = subparsers.add_parser('JCOUPLING', help='Calculate J couplings.')
    parser_09.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_09.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_09.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_09.add_argument("--top", dest="top",help="Topology file",required=False)
    parser_09.add_argument("--res", dest="res",help="Calculate couplings for specific residues",required=False,nargs="+",default=None)
    parser_09.add_argument("--raw",dest="raw",help="print raw angles for j3",action="store_true",default=False)

    # CREATE FRAGMENTS
    parser_07 = subparsers.add_parser('SNIPPET', help='SPLIT structure in multiple PDB')
    parser_07.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_07.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=True)
    parser_07.add_argument("--seq", dest="seq",help="Sequence type. Accepts ACGU/NRY format",required=True)
    parser_07.add_argument("--outdir", dest="outdir",help="outdir",required=False,default=None)

    # CALCULATE ENM
    parser_10 = subparsers.add_parser('ENM', help='Calculate ENM')
    parser_10.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_10.add_argument("--pdb", dest="pdbs",help="PDB file",default=None,required=True)

    parser_10.add_argument("--cutoff", dest="cutoff",help="Cutoff distance in nm (default=0.9)",default=0.9,type=float)
    parser_10.add_argument("--type", dest="type",default='SBP',choices=['P','S','B','SB','SP','BP','SBP','AA'], help='Type of ENM (default=SBP)')    
    parser_10.add_argument("--ntop", dest="ntop",help="Number of top eigenvectors to write (default=10)",default=10,type=int)
    parser_10.add_argument("--sparse", dest="sparse",help="Use sparse matrix diagonalization algorithm. Recomended for large matrices.",action='store_true',default=False)

    parser_10.add_argument("--shape", dest="shape",help="calculate C2/C2 fluctuations",action='store_true',default=False)
    parser_10.add_argument("--modes", dest="modes",help="dump top modes on a trajectory",action='store_true',default=False)

    
    # CLUSTERING
    #parser_13 = subparsers.add_parser('CLUSTER', help="Clustering")
    #parser_13.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    #parser_13.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    #parser_13.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    #parser_13.add_argument("--top", dest="top",help="Topology file",required=False)
    #parser_13.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    
    #parser_13.add_argument("--alg", dest="alg",required=False,default="DBSCAN",choices=["DBSCAN"])
    #parser_13.add_argument("--eps", dest="eps",help="epsilon used for DBSCAN clustering. default=0.7",default=0.7,type=float)
    #parser_13.add_argument("--mins", dest="mins",help="minsample for DBSCAN clustering, Default=50",default=50,type=float)
    
    args = parser.parse_args(arguments)

    # CLUSTER
    
    return args


####################### ERMSD ########################

def ermsd(args):

    if(args.top==None):
        dd = [bb.ermsd(args.reference,pdb,cutoff=args.cutoff) for pdb in args.pdbs]
    else:
        dd = bb.ermsd(args.reference,args.trj,topology=args.top,cutoff=args.cutoff)
        
    fh = open(args.name + ".out",'w')
    fh.write("# %s \n" % (" ".join(sys.argv[:])))
    fh.write("#%10s   %10s\n" % ("Frame","eRMSD"))
    fh.write("".join([ " %10d   %10.4e \n" % (i,d) for i,d in enumerate(dd)]))
    fh.close()

####################### RMSD ########################

def rmsd(args):


    if(args.top==None):
        if(args.dump==True):
            dd = [bb.rmsd(args.reference,pdb,out="%s_%06d.pdb" % (args.name,i),heavy_atom=args.heavy_atom) for i,pdb in enumerate(args.pdbs)]
        else:
            dd = [bb.rmsd(args.reference,pdb,heavy_atom=args.heavy_atom) for i,pdb in enumerate(args.pdbs)]
    else:
        if(args.dump==True):
            out = "%s.%s" % (args.name, (args.trj).split(".")[-1])
            dd = bb.rmsd(args.reference,args.trj,topology=args.top,out=out,heavy_atom=args.heavy_atom)
        else:
            dd = bb.rmsd(args.reference,args.trj,topology=args.top,heavy_atom=args.heavy_atom)

    fh = open(args.name + ".out",'w')
    fh.write("# %s \n" % (" ".join(sys.argv[:])))
    fh.write("#%10s   %10s\n" % ("Frame","RMSD"))
    fh.write("".join([ " %10d   %10.4e \n" % (i,d) for i,d in enumerate(dd)]))
    fh.close()

####################### ESCORE ########################
        
def score(args):
    
    import barnaba.escore as escore
    # set "force-field"
    ee = escore.Escore([args.reference],cutoff=args.cutoff)
    if(args.top==None):
        dd = [ee.score(pdb)[0] for pdb in args.pdbs]
    else:
        dd = ee.score(args.trj,topology=args.top)

    # Write to file
    fh = open(args.name + ".out",'w')
    fh.write("# %s \n" % (" ".join(sys.argv[:])))
    fh.write("#%10s   %10s\n" % ("Frame","ESCORE"))
    fh.write("".join([ " %10d   %10.4e \n" % (i,d) for i,d in enumerate(dd)]))
    fh.close()


####################### SS_MOTIF ########################

def ss_motif(args):

    out = None
    if(args.dump==True):out = args.name

    stri = "# %s \n" % (" ".join(sys.argv[:]))
    if(args.top==None):
        stri += "#%-20s %10s %s \n" % ("PDB","eRMSD","Sequence")
        for i in range(len(args.pdbs)):
            try:
                dd = bb.ss_motif(args.query,args.pdbs[i],out=out,bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
                stri += "".join([" %-20s %10.4e %s \n" % (args.pdbs[i].split("/")[-1],dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])
            except:
                print("# not able to load %s" % args.pdbs[i])
                continue
    else:
        stri += "#%-10s %-10s %10s %s \n" % ("index","frame","eRMSD","Sequence")
        dd = bb.ss_motif(args.query,args.trj,topology=args.top,out=out,bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
        stri += "".join([" %-10d %-10d %10.4e %s \n" % (j,dd[j][0],dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])

    fh = open(args.name + ".out",'w')
    fh.write(stri)
    fh.close()


####################### DS_MOTIF ########################

def ds_motif(args):

    out = None
    if(args.dump==True):out = args.name

    stri = "# %s \n" % (" ".join(sys.argv[:]))
    if(args.top==None):
        stri += "#%-20s %10s %s \n" % ("PDB","eRMSD","Sequence")
        for i in range(len(args.pdbs)):
            dd = bb.ds_motif(args.query,args.pdbs[i],out=out,l1=args.l1,l2=args.l2,\
                             bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
            stri += "".join([" %-20s %10.4e %s \n" % (args.pdbs[i].split("/")[-1],dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])
    else:
        stri += "#%-10s %-10s %10s %s \n" % ("index","frame","eRMSD","Sequence")
        dd = bb.ss_motif(args.query,args.trj,topology=args.top,out=out,\
                         bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
        #stri += "".join([" %-20d %10.4e %s \n" % (j,dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])
        stri += "".join([" %-10d %-10d %10.4e %s \n" % (j,dd[j][0],dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])

    fh = open(args.name + ".out",'w')
    fh.write(stri)
    fh.close()


####################### ANNOTATE ########################

def annotate(args):


    stri_p = "# %s \n" % (" ".join(sys.argv[:]))
    stri_p += "#%-10s %-10s %4s \n" % ("RES1","RES2","ANNO")
    
    stri_s = "# %s \n" % (" ".join(sys.argv[:]))
    stri_s += "#%-10s %-10s %4s \n" % ("RES1","RES2","ANNO")

    stri_dot = "# %s \n" % (" ".join(sys.argv[:]))
    
    if(args.top==None):
        for i in range(len(args.pdbs)):
            st, pair, res = bb.annotate(args.pdbs[i])
            stri_p += "# PDB %s \n" % args.pdbs[i].split("/")[-1]
            stri_s += "# sequence %s\n" % "-".join(res)
            stri_p += "# sequence %s\n" % "-".join(res)
            stri_p += "".join([ "%-10s %-10s %4s \n" % (res[pair[0][0][e][0]],res[pair[0][0][e][1]],pair[0][1][e]) for e in range(len(pair[0][0]))])
            
            stri_s += "# PDB %s \n" % args.pdbs[i].split("/")[-1]
            stri_s += "".join([ "%-10s %-10s %4s \n" % (res[st[0][0][e][0]],res[st[0][0][e][1]],st[0][1][e]) for e in range(len(st[0][0]))])

            if(args.dotbr):
                dotbr,ss = bb.dot_bracket(pair,res)
                stri_dot += "# sequence %s\n" % "-".join(res)
                stri_dot += "# SEQ %s" % ss
                stri_dot += "%-20s %s\n" %(args.pdbs[i].split("/")[-1],dotbr[0])

            
    else:
        st,pair,res = bb.annotate(args.trj,topology=args.top)
        stri_s += "# sequence %s\n" % "-".join(res)
        stri_p += "# sequence %s\n" % "-".join(res)
        stri_dot += "# sequence %s\n" % "-".join(res)

        if(args.dotbr):
            dotbr,ss = bb.dot_bracket(pair,res)
            stri_dot += "# SEQ %s" % ss
            stri_dot += " ".join(["%-10d %s\n" %(k,dotbr[k]) for k in range(len(pair))])
            
        for k in range(len(st)):
            stri_p += "# Frame %d \n" % k
            stri_p += "".join([ "%-10s %-10s %4s \n" % (res[pair[k][0][e][0]],res[pair[k][0][e][1]],pair[k][1][e]) for e in range(len(pair[k][0]))])
        
            stri_s += "# Frame %d \n" % k
            stri_s += "".join([ "%-10s %-10s %4s \n" % (res[st[k][0][e][0]],res[st[k][0][e][1]],st[k][1][e]) for e in range(len(st[k][0]))])
            
            
    fh1 = open(args.name + ".pairing.out",'w')
    fh1.write(stri_p)
    fh1.close()
    
    fh2 = open(args.name + ".stacking.out",'w')
    fh2.write(stri_s)
    fh2.close()

    if(args.dotbr):
        fh3 = open(args.name + ".dotbracket.out",'w')
        fh3.write(stri_dot)
        fh3.close()


####################### SEC_STRUCTURE #####################

def sec_structure(args):
    args.chi_conf = 0
    import datetime
    start_time = datetime.datetime.now()
    import numpy as np
  #  from numpy import unravel_index
    import barnaba.sec_str_svg as sesvg
    import barnaba.sec_str_ff as seff
    import barnaba.sec_str_constants as secon

    output = "# %s \n" % (" ".join(sys.argv[:]))
    weights = []
    if args.not_order:
        order = False
    else:
        order = True
    threshold = float(args.thresh)
#    if args.f_weights:
#        with open(args.f_weights) as f:
#            content = f.readlines()
#            f.close()
#        for line in content:
#            if not line.startswith("#"):
#                weights.append(float(line.split()[-1]))
#        print("Found weights file %s with %d weights" % (args.f_weights, n_weights))
    n_weights = len(weights)
    


    a_ann_list = {}
    a_pairs = []
    n_frames = 0
    f_inp = {}
    for f in args.f_anns:
        if f.endswith("ANNOTATE.dotbracket.out"):
            f_inp['db'] = f
        elif f.endswith("ANNOTATE.pairing.out"):
            f_inp['pa'] = f
        elif f.endswith("ANNOTATE.stacking.out"):
            f_inp['st'] = f
        else:
            print("Cannot recognize annotation file, try pairing.")
            f_inp['pa'] = f

    if 'db' in f_inp :
        if 'pa' in f_inp or 'st' in f_inp:
            if not args.template:
                print("Using dotbracket for structure minimization, but showing all annotations.")
            else:
                print("Using template structure with barnaba annotations.")
            mode = 2
        else:
            if not args.template:
                print("Using dotbracket for structure minimization.")
            else:
                print("Using template structure with barnaba dotbracket.")
            mode = 1
    elif 'pa' in f_inp or 'st' in f_inp:
            if not args.template:
                print("Using annotation file(s) for structure minimization.")
            else:
                print("Using template structure with barnaba annotations.")
            mode = 0

    for k, f in f_inp.items():
        if k == 'db':
            db_sequence, db_ann_list, db_pairs, db_n_frames = bb.parse_dotbracket(threshold, f, weights)
       #     pairs_list.append(db_pairs)
            dotbr = True
        else:    
            i_sequence, i_ann_list, i_pairs, i_n_frames = bb.parse_annotations(threshold, f, weights)
            if len(a_ann_list) == 0: a_ann_list = i_ann_list
            else: a_ann_list.update(i_ann_list)
            if len(a_pairs) == 0: a_pairs = i_pairs
            else: a_pairs = np.unique(np.append(a_pairs, i_pairs, axis=0), axis=0)
            if n_frames == 0:
                n_frames = i_n_frames
            else:
                if n_frames != i_n_frames:
                    sys.exit("Annotations files %s have different numbers of frames" % args.f_anns)
                if n_weights > 0 and n_weights != n_frames:
                    sys.exit("Annotations files %s have different numbers of frames than number of weights" % args.f_anns)
    if mode == 2 and not (db_sequence==i_sequence).all():
        sys.exit("Sequences in dotbracket and anntation file(s) differ.")
    if mode == 2 and not (db_n_frames == i_n_frames):
        sys.exit("Number of frames in dotbracket and anntation file(s) differ.")
        

#    if len(chains) > 1:
#        output += "# Found >1 chains, but can draw only single strand.\n# Drawing first chain.\n"
#        print("# Found >1 chains, but can draw only single strand.\n# Drawing first chain.")

    c = 0
    if mode in [1,2]:
        ann_list = db_ann_list 
        pairs = np.array(db_pairs)
        sequence = db_sequence
    else:
        ann_list = a_ann_list 
        pairs = a_pairs
        sequence = i_sequence
    if len(ann_list) == 0:
        sys.exit("No annotations found.")

    if args.template:
        import re
        with open(args.template) as f:
            content = f.readlines()
            f.close()
        re_position_varna = re.compile("<circle cx=\"([0-9\.]+)\" cy=\"([0-9\.]+)\" r=\"5.0\" stroke=\"none\" stroke-width=\"1.0\" fill=\"rgb\(.*\)\"/>")
        re_position_rnastructure = re.compile(".*<circle style=\"stroke-width:1\" cx=\"([0-9\.]+)\" cy=\"([0-9\.]+)\" .*/>")
        ri = 0
        pos = np.zeros((len(sequence), 2))
        for re_position in [re_position_varna, re_position_rnastructure]:
            for line in content:
                _pos = re_position.match(line)
                if _pos:
                    x = float(_pos.group(1))
                    y = float(_pos.group(2))
                    pos[ri] = [x, y]
                    ri += 1        
            if ri == 0:
                continue
            if ri != len(sequence):
                print(ri, len(sequence))
                sys.exit("Length of annotation sequence is different of number of residues in template.")
            else:
                break
        scale = secon.d_seq/np.mean(np.linalg.norm(pos[1:]-pos[:-1], axis=1))
        pos *= scale
        dimensions = np.linalg.norm([pos[:,0].max()-pos[:,0].min(), pos[:,1].max()-pos[:,1].min()])
        if mode == 2:
            pairs = a_pairs
            ann_list = a_ann_list
        output_svg = sesvg.draw_structure(threshold, pos, pairs, ann_list, int(args.chi_conf), sequence, dimensions, args.output_ids)
        fh1 = open(args.name + "_from_template.svg",'w')
        fh1.write(output_svg)
        fh1.close()
        sys.exit()
        
    n = len(sequence)
    __i1, __i2 = np.triu_indices(n, k=1)
    __k_rep = np.full((int(n*(n-1)/2), 2), secon.k_rep2)
    __k_rep_lr = np.full((int(n*(n-1)/2), 2), secon.k_rep_lr)
    __d_rep = np.full((int(n*(n-1)/2), 2), secon.d_rep2)
    print("Calc parameters ...")
    param_seq, param_rep, param_rep_lr, sorted_params, param_stem, param_ang, param_bulge, param_bulge_rep, param_angle_180 = bb.parameters(pairs, ann_list, n, threshold)
    print("... done")
 #   print(param_seq)
 #   print(param_rep)
 #   print(sorted_params)
 #   print(param_stem)
 #   print(param_ang)
 #   print(param_angle_180)
 #   sys.exit()
    if len(sorted_params) == 0:
        sys.exit("No annotations found. Exiting.")        
    dimensions = secon.d_seq*(n-1)*.7
    length = n*secon.d_seq
    angle = 2*np.pi/n
    start = np.column_stack((np.arange(0., n*secon.d_seq, secon.d_seq),  np.full((n), 0.5*dimensions)))
    i_min = n-1
    i_max = 0

    term_pair = sorted_params[np.argmax(sorted_params[:,2]-sorted_params[:,1])]
    cotranscript = False
    param = {}
    for p in sorted_params:
        if ((min(p[1:3]) <  min(term_pair[1:3]) < max(p[1:3])) or
             ((min(p[1:3]) <  max(term_pair[1:3]) < max(p[1:3])))):
            cotranscript = True
            break
#    cotranscript = True        
    param[0] = param_seq
    param[1] = param_rep
    steps = int(args.nsteps)
    if not cotranscript:
        param[0] = np.append(param[0], sorted_params, axis=0) 
        print("Start on circle")
    print("Make starting config")
    if cotranscript:
       # start on straight slightly perturbed line
       p, ii_pair = bb.get_par(sorted_params, start)
       for i in np.arange(p[1]+1, p[2]):
            start[int(i), 1] -= 5 
    else:    
        for i in range(n):
       # start as circle
            start[i][0] = -length * 0.5/np.pi * np.sin(i*angle + .5*angle) + dimensions * 0.5
            start[i][1] = length * 0.5/np.pi * np.cos(i*angle + .5*angle) + dimensions * 0.5
       # start on half circle
       # start[i][0] = length * 1./np.pi * np.sin(i*angle*.5 + .5*angle) + dimensions * 0.5
       # start[i][1] = length * 1./np.pi * np.cos(i*angle*.5 + .5*angle) + dimensions * 0.5
    pos = start
    h = secon.h
    print_snapshots = min([steps, int(args.draw_interm)])
    print_energy = 20    
    print_energy = min([steps, print_energy])
    ds = int(float(steps)/print_energy)
    if print_snapshots > 0:
        ds_draw = int(float(steps)/print_snapshots)
    else:
        ds_draw = 0    
#       dt = float(args.T_init)/int(steps)*3/2
#       T = float(args.T_init)
        
    write_force = False
    i_pair = 0
    i_ang =  0
    i_par = 0
    param[7] = param_rep_lr
    param[8] = param_angle_180
    E = seff.energy(pos, param)
    new_E = E
    F = seff.force(pos, param, __i1, __i2, __k_rep, __d_rep, __k_rep_lr, write_force)
    max_F = np.linalg.norm(F, axis=1).max()
    r_i = np.linalg.norm(F, axis=1).argmax()
    E_thresh = 5e-5
    output += "# E thresh: %.2e\n" % E_thresh
    print("# E thresh: %.2e" % E_thresh)
    output += "%8s%10s%10s%10s%10s\n" % ("Step", "energy", "F_max",  "h", "res(max_F)") 
    print("%8s%10s%10s%10s%10s" % ("Step", "energy", "F_max",  "h", "res(max_F)"))
    if order:
        reduce_straight_pot = True
    else:    
        reduce_straight_pot = False
#    print("Prep. took", datetime.datetime.now()-start_time)
    sorted_params_copy = np.copy(sorted_params)
    added_ang_pars = []
    added_sort_params = np.empty((0, 5))
    for i in range(100000):
        print_out = False
        if (new_E <= E and  (E-new_E)/E < E_thresh) or h<1e-7:
            if cotranscript and len(sorted_params) > 0:    
                p, ii_pair = bb.get_par(sorted_params, pos)
                print_out = True
                param[0] = np.append(param[0], [p], axis=0)
                output += "%d Added %s\n" % (i, p)
                print(i, " Added ", p, len(sorted_params))
                sorted_params = np.delete(sorted_params, ii_pair, 0)
                if len(added_sort_params) == 0:
                    added_sort_params = np.array([p])
                else:
                    added_sort_params = np.append(added_sort_params, [p], axis=0)
                if len(added_ang_pars) < len(param_ang):
                    keys1, d = np.where(param_ang[:,1:3]==p[1:3])
                    keys2, d = np.where(param_ang[:,3:5]==p[1:3])
                    if len(keys1) == 2 and keys1[0] not in added_ang_pars:
                        i_ang = keys1[0]
                        if param_ang[i_ang, 3:5] in added_sort_params[:, 1:3]:
                            if len(added_ang_pars) == 0:
                                param[2] = np.array([param_ang[i_ang]])
                            else: param[2] = np.append(param[2], [param_ang[i_ang]], axis=0)
                            print("Added soft ang. pot.:", param_ang[i_ang,1:5])
                            added_ang_pars.append(i_ang)
                    if len(keys2) == 2 and keys2[0] not in added_ang_pars:
                        i_ang = keys2[0]
                        if param_ang[i_ang, 1:3] in added_sort_params[:, 1:3]:
                            if len(added_ang_pars) == 0:
                                param[2] = np.array([param_ang[i_ang]])
                            else: param[2] = np.append(param[2], [param_ang[i_ang]], axis=0)
                            added_ang_pars.append(i_ang)
                            print("Added soft ang. pot.:", param_ang[i_ang,1:5])
                i_ang = 0      
                E = seff.energy(pos, param)
                h = secon.h
                if len(sorted_params) == 0 and steps-i < 500:
                    steps = i+500
                    output += "Setting new step limit %d\n" % steps
                    print("Setting new step limit %d" % steps)
            elif order:
                print_out = True
                if i_ang == 0:
                    param_ang[:,5] *= secon.k_ang_end/secon.k_ang
                    if not cotranscript: param[2] = np.array([param_ang[i_ang]])
                    else: param[2] = np.append(param[2], [param_ang[i_ang]], axis=0)
                else:    
                    param[2] = np.append(param[2], [param_ang[i_ang]], axis=0)
                output += "%d Added ordering pot %s\n" % (i, param_ang[i_ang,1:5])    
                print("Added ordering pot", param_ang[i_ang,1:5])
                h = secon.h
                i_ang += 1
                E = seff.energy(pos, param)
                if i_ang == len(param_ang): 
                    order = False
                    
            elif not order and reduce_straight_pot:
                __k_rep_lr = np.full((int(n*(n-1)/2), 2), secon.k_rep_lr_end)
                print_out = True
                p_tmp = param_angle_180
                p_tmp[:,4] = secon.k_angle_straight_end
                param[8] = p_tmp
                E = seff.energy(pos, param)
                output += "%d Relax straightening potential\n" % i
                print("Reduced straight pot")
                reduce_straight_pot = False
                if steps-i < 500:
                    steps = i+500
                    output += "Setting new step limit %d\n" % steps
                    print("Setting new step limit %d" % steps)
                h = secon.h
            else:
                E = new_E
        else:
            E = new_E

        if i % ds == 0 or print_out:
            r_i = np.linalg.norm(F, axis=1).argmax()
            print(    "%8d%10.3e%10.2e%10.2e%10d" % (i, E, max_F, h,  r_i))
            output += "%8d%10.2e%10.2e%10.2e%10d\n" % (i, E, max_F, h, r_i) 
            write_force = False

        else:
            write_force = False
        if ds_draw > 0 and i % ds_draw == 0:
            if mode in [0, 2]:
                ann_list = a_ann_list    
                pairs = a_pairs
            output_svg = sesvg.draw_structure(threshold, pos, pairs, ann_list, int(args.chi_conf), sequence, dimensions, args.output_ids)
            fh1 = open(args.name + "_%03d.svg" % i,'w')
            fh1.write(output_svg)
            fh1.close()
            
        if i == steps:
            output += "Reached maximum number of %d steps\n" % i
            print("Reached maximum number of %d steps\n" % i)
            if order or reduce_straight_pot:
                print("Ordering potential has not yet been added, ignoring step limit")
                output += "Ordering potential has not yet been added, ignoring step limit\n"
            elif (cotranscript and len(sorted_params) > 0):    
                print("Not all interactions have been added, ignoring step limit")
                output += "Not all interactions have been added, ignoring step limit\n"
            else:    
                break
        F = seff.force(pos, param, __i1, __i2, __k_rep, __d_rep, __k_rep_lr, write_force)
        max_F = np.linalg.norm(F, axis=1).max()
        new_pos = pos + F/max_F * h
        new_E = seff.energy(new_pos, param)
 
        while new_E >= E and h > 1e-8 and ((order or reduce_straight_pot) or h > 1e-4):
            h *= 0.2
            new_pos = pos + F/max_F * h
            new_E = seff.energy(new_pos, param)
        h *= 1.2

        if h < 1e-4 and not order and not reduce_straight_pot and not (cotranscript and len(sorted_params) > 0): 
            output += "Converged at %d steps.\n" % i
            print("Converged at %d steps.\n" % i)
            break
        pos = new_pos-np.full((len(pos), 2), (np.mean(new_pos, axis=0)-[dimensions*.5, dimensions*.5]))
        if h > secon.h:
            h = secon.h
            
    E = seff.energy(pos, param, False) 
    F = seff.force(pos, param, __i1, __i2, __k_rep, __d_rep, __k_rep_lr, write_force)
    r_i = np.linalg.norm(F, axis=1).argmax()
    print("%8d%10.2e%10.2e%10.2e%10d" % (i, E, max_F,  h, r_i))
    output += "%8d%10.2e%10.2e%10.2e%10d\n" % (i, E, max_F, h, r_i) 
    print( i, " steps of minmization")    
    output += "%d steps of minmization\n" % (i)
    from scipy.linalg import expm, norm

    av_stem_vector = np.mean([pos[int(param_stem[0][1])], pos[int(param_stem[0][2])]], axis=0) - np.mean([pos[int(param_stem[len(param_stem)-1][1])], pos[int(param_stem[len(param_stem)-1][2])]], axis=0)
    theta = np.arctan2(-av_stem_vector[0], -av_stem_vector[1])

    def M(axis, theta):
        return expm(np.cross(np.eye(3), axis/norm(axis)*theta))
    axis = [0,0,1]
    m0 = M(axis, theta)
    trans_pos = np.full((len(pos),3), [dimensions*.5, dimensions*.5, 0])
    pos = np.hstack((pos,np.zeros((pos.shape[0],1))))
    pos  = np.array(((np.asmatrix(m0) * np.asmatrix(pos-trans_pos).T).T+trans_pos)[:,0:2])
    if pos[0][0] > pos[len(sequence)-1][0]:
        m0 = M(axis, np.pi)
        trans_pos = np.full((len(pos),3), [dimensions*.5, dimensions*.5, 0])
        pos = np.hstack((pos,np.zeros((pos.shape[0],1))))
        pos  = np.array(((np.asmatrix(m0) * np.asmatrix(pos-trans_pos).T).T+trans_pos)[:,0:2])
    if mode in [0, 2]:
        ann_list = a_ann_list
        pairs = a_pairs
    output_svg = sesvg.draw_structure(threshold, pos, pairs, ann_list, int(args.chi_conf), sequence, dimensions, args.output_ids)
    fh1 = open(args.name + "_%dsteps.svg" % (i),'w')
    fh1.write(output_svg)
    fh1.close()


    fh2 = open(args.name + ".out",'w')
    fh2.write(output)
    end_time = datetime.datetime.now()
    fh2.write("# Took %s\n" % (end_time - start_time))
    fh2.close()

    
                     



##################### DUMP #######################

def dump(args):

    assert args.dumpR or args.dumpG, "# ERROR. choose --dumpR and/or --dumpG"
    

    if(args.dumpR):
        stri_r = "# %s \n" % (" ".join(sys.argv[:]))
        stri_r += "#%15s %15s %11s %11s %11s \n" % ("RES1","RES2","x","y","z")

        if(args.top==None):
            for i in range(len(args.pdbs)):
                rvecs,resi = bb.dump_rvec(args.pdbs[i],cutoff=args.cutoff)
                idxs = its.permutations(range(len(resi)), 2)
                stri_r += "# PDB %s \n" % args.pdbs[i].split("/")[-1]
                stri_r += "".join([" %15s %15s %11.4e %11.4e %11.4e \n" % (resi[i1],resi[i2],rvecs[0,i1,i2,0],rvecs[0,i1,i2,1],rvecs[0,i1,i2,2]) for i1,i2 in idxs if(sum(rvecs[0,i1,i2]**2)> 1.E-05)])
        else:
            rvecs,resi = bb.dump_rvec(args.trj,topology=args.top,cutoff=args.cutoff)
            idxs = its.permutations(range(len(resi)), 2)
            for i in range(len(rvecs)):
                stri_r += "# Frame %d \n" % i
                stri_r += "".join([" %15s %15s %11.4e %11.4e %11.4e \n" % (resi[i1],resi[i2],rvecs[i,i1,i2,0],rvecs[i,i1,i2,1],rvecs[i,i1,i2,2]) for i1,i2 in idxs if(sum(rvecs[i,i1,i2]**2)> 1.E-05)])
                
        fh = open(args.name + ".rvec.out",'w')
        fh.write(stri_r)
        fh.close()
        
    if(args.dumpG):
        
        stri_g = "# %s \n" % (" ".join(sys.argv[:]))
        stri_g += "#%15s %15s %11s %11s %11s %11s \n" % ("RES1","RES2","G0","G1","G2","G3")
        
        if(args.top==None):
            for i in range(len(args.pdbs)):
                rvecs,resi = bb.dump_gvec(args.pdbs[i],cutoff=args.cutoff)
                idxs = its.permutations(range(len(resi)), 2)
                stri_g += "# PDB %s \n" % args.pdbs[i].split("/")[-1]
                stri_g += "".join([" %15s %15s %11.4e %11.4e %11.4e %11.4e \n" % (resi[i1],resi[i2],rvecs[0,i1,i2,0],rvecs[0,i1,i2,1],rvecs[0,i1,i2,2],rvecs[0,i1,i2,3]) for i1,i2 in idxs if(sum(rvecs[0,i1,i2]**2)> 1.E-05)])
        else:
            rvecs,resi = bb.dump_rvec(args.trj,topology=args.top,cutoff=args.cutoff)
            idxs = its.permutations(range(len(resi)), 2)
            for i in range(len(rvecs)):
                stri_g += "# Frame %d \n" % i
                stri_g += "".join([" %15s %15s %11.4e %11.4e %11.4e %11.4e \n" % (resi[i1],resi[i2],rvecs[i,i1,i2,0],rvecs[i,i1,i2,1],rvecs[i,i1,i2,2],rvecs[i,i1,i2,3]) for i1,i2 in idxs if(sum(rvecs[i,i1,i2]**2)> 1.E-05)])
                
        fh = open(args.name + ".gvec.out",'w')
        fh.write(stri_g)
        fh.close()
        

##################### CALCULATE TORSION ANGLES #######################
    
def torsion(args):
    
    assert args.backbone or args.sugar or args.pucker, "# ERROR. choose --backbone/sugar/pucker"
    
    if(args.backbone):
        
        stri_b = "# %s \n" % (" ".join(sys.argv[:]))        
        stri_b += "#%-12s  %11s %11s %11s %11s %11s %11s %11s\n" % ("RESIDUE","alpha","beta","gamma","delta","eps","zeta","chi")

        if(args.top==None):
            for i in range(len(args.pdbs)):
                stri_b += "# PDB %s \n" % args.pdbs[i].split("/")[-1]                
                angles_b,rr = bb.backbone_angles(args.pdbs[i],residues=args.res)
                stri_b += "".join([" %-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[0,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
        else:
            
            angles_b,rr = bb.backbone_angles(args.trj,topology=args.top,residues=args.res)
            for i in range(angles_b.shape[0]):
                stri_b += "# Frame %d \n" % i
                stri_b += "".join([" %-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[i,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
            
            
        fh = open(args.name + ".backbone.out",'w')
        fh.write(stri_b)
        fh.close()

    if(args.sugar):
        
        stri_b = "# %s \n" % (" ".join(sys.argv[:]))        
        stri_b += "#%-12s  %11s %11s %11s %11s %11s \n" % ("RESIDUE","nu0","nu1","nu2","nu3","nu4")

        if(args.top==None):
            for i in range(len(args.pdbs)):
                stri_b += "# PDB %s \n" % args.pdbs[i].split("/")[-1]                
                angles_b,rr = bb.sugar_angles(args.pdbs[i],residues=args.res)
                stri_b += "".join([" %-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[0,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
        else:
            
            angles_b,rr = bb.sugar_angles(args.trj,topology=args.top,residues=args.res)
            for i in range(angles_b.shape[0]):
                stri_b += "# Frame %d \n" % i
                stri_b += "".join([" %-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[i,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
            
            
        fh = open(args.name + ".sugar.out",'w')
        fh.write(stri_b)
        fh.close()

    if(args.pucker):
        
        stri_b = "# %s \n" % (" ".join(sys.argv[:]))        
        stri_b += "#%-12s  %11s %11s \n" % ("RESIDUE","Phase","Amplitude")

        if(args.top==None):
            for i in range(len(args.pdbs)):
                stri_b += "# PDB %s \n" % args.pdbs[i].split("/")[-1]                
                angles_b,rr = bb.pucker_angles(args.pdbs[i],residues=args.res,altona=args.altona)
                stri_b += "".join(["%-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[0,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
        else:
            
            angles_b,rr = bb.pucker_angles(args.trj,topology=args.top,residues=args.res,altona=args.altona)
            for i in range(angles_b.shape[0]):
                stri_b += "# Frame %d \n" % i
                stri_b += "".join(["%-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[i,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
                        
        fh = open(args.name + ".pucker.out",'w')
        fh.write(stri_b)
        fh.close()
      

def couplings(args):

    from  barnaba import definitions
    stri = "# %s \n" % (" ".join(sys.argv[:]))
    stri += "#%-12s %s\n" % ("RESIDUE","".join([" %11s" % (k) for k in  definitions.couplings_idx.keys()]))
    
    if(args.top==None):
        for i in range(len(args.pdbs)):
            stri += "# PDB %s \n" % args.pdbs[i].split("/")[-1]                
            angles_b,rr = bb.jcouplings(args.pdbs[i],residues=args.res,raw=args.raw)
            stri += "".join(["%-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[0,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
    else:
        angles_b,rr = bb.jcouplings(args.trj,topology=args.top,residues=args.res,raw=args.raw)
        for i in range(angles_b.shape[0]):
            stri += "# Frame %d \n" % i
            stri += "".join([" %-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[i,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])

    fh = open(args.name + ".couplings.out",'w')
    fh.write(stri)
    fh.close()

    

####################  SNIPPET #######################

def snippet(args):

    for pdb in args.pdbs:
        bb.snippet(pdb,args.seq,outdir=args.outdir)


####################  ENM  #######################

def enm(args):
    
    import barnaba.enm as enm
    
    if(args.type=="AA"):
        sele = "AA"
    else:
        sele = []
        sele.append("CA")
        if("P" in args.type):
            sele.append("P")
        if("S" in args.type):
            sele.append("C1\'")
        if("B" in args.type):
            sele.append("C2")
        if(args.type=="SBP"):
            sele.append("CB")

    net = enm.Enm(args.pdbs,sele_atoms=sele,sparse=args.sparse,ntop=args.ntop,cutoff=args.cutoff)

    # print eigenvectors 
    eigvecs = net.print_evec(args.ntop)
    fh = open(args.name + ".eigvecs.out",'w')
    fh.write(eigvecs)
    fh.close()

    # print eigenvalues 
    eigvals = net.print_eval()
    fh = open(args.name + ".eigvals.out",'w')
    fh.write(eigvals)
    fh.close()

    if(args.shape):
        
        fluc,res =  net.c2_fluctuations()
        stri = "# %19s %s \n" % ("Residues","Fluctuations")
        stri +=  "".join(["%10s/%-10s %.6e \n" % (res[i],res[i+1],fluc[i]) for i in range(len(fluc)) ])
        fh = open(args.name + ".shape.out",'w')
        fh.write(stri)
        fh.close()

    if(args.modes):
        for i in range(6,args.ntop+6):
           mode=net.get_mode_traj(i)
           mode.save_pdb(args.name + ".mode_" + str(i) + ".pdb")


####################### MAIN #########################

def main(arguments=None):
    """
    Use barnaba command line tool

    Parameters
    ----------
    arguments :
         List of command line arguments. If not provided, positional arguments will be used.

    Returns error code. 0 means no error.
    """

# This is to allow passing a single string
    if(isinstance(arguments,str)):
        arguments=arguments.split()

    def filename(args):
        # create output filename
        if(args.name == None):
            outfile = 'outfile.' + args.subparser_name 
        else:
            outfile = args.name + "." + args.subparser_name
        args.name = outfile
        print("# your output will be written to files with prefix %s" % outfile)

    # Parse options
    # In order to avoid python to crash when wrong arguments are used, it is
    # necessary to intercept the exception.
    try:
        args = parse(arguments)
    except:
        return sys.exc_info()[0]

    # create filename
    if(args.subparser_name!="SNIPPET"):
        filename(args)

    # check
    if(args.subparser_name!="SEC_STRUCTURE"):
        if(args.pdbs==None):
            assert args.trj != None, "# Specify either pdbs (--pdb) or a trajectory file (--trj)"
            assert args.top != None, "# Please provide a topology file"
        else:
            if(args.subparser_name != "ENM" and args.subparser_name!="SNIPPET"):
                assert args.trj == None, "# Specify either pdbs (--pdb) or a trajectory+topology files, not both"
            
    
    # call appropriate function
    options = {'ERMSD' : ermsd,\
               'RMSD' : rmsd,\
               'ESCORE' : score,\
               'SS_MOTIF' : ss_motif,\
               'DS_MOTIF' : ds_motif,\
               'ANNOTATE' : annotate,\
               'SEC_STRUCTURE' : sec_structure,\
               'DUMP' : dump,\
               'TORSION':torsion,\
               'JCOUPLING':couplings,\
               'ENM':enm,\
               'SNIPPET' : snippet}

    options[args.subparser_name](args)
    
