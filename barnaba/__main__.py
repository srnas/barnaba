#!/usr/bin/env python

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

def parse():

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

    # SEARCH DOUBLE STRANDED MOTIFS - OK
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

    parser_08.add_argument("--backbone", dest="backbone",help="calculate backbone (a,b,g,d,e,z) and chi torsion angle",action='store_true',default=False)
    parser_08.add_argument("--sugar", dest="sugar",help="calculate sugar torsion angles (v0...v5)",action='store_true',default=False)
    parser_08.add_argument("--pucker", dest="pucker",help="calculate sugar pucker angles (amplitude, phase)",action='store_true',default=False)
    parser_08.add_argument("--res", dest="res",help="Calculate torsion angle for specific residues",required=False,nargs="+",default=None)
    

    # CALCULATE TORSION ANGLES - J COUPLINGS OK
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
    
    args = parser.parse_args()

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
            dd = [bb.rmsd(args.reference,pdb,out="%s_%06d.pdb" % (args.name,i)) for i,pdb in enumerate(args.pdbs)]
        else:
            dd = [bb.rmsd(args.reference,pdb) for i,pdb in enumerate(args.pdbs)]
    else:
        if(args.dump==True):
            out = "%s.%s" % (args.name, (args.trj).split(".")[-1])
            dd = bb.rmsd(args.reference,args.trj,topology=args.top,out=out)
        else:
            dd = bb.rmsd(args.reference,args.trj,topology=args.top)

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
        dd = bb.ss_motif(args.query,args.trj,topology=args.top,out=out,bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
        stri += "".join([" %-20d %10.4e %s \n" % (j,dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])

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
        stri += "#%20s %10s %s \n" % ("frame","eRMSD","Sequence")
        dd = bb.ss_motif(args.query,args.trj,topology=args.top,out=out,l1=args.l1,l2=args.l2,\
                         bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
        stri += "".join([" %-20d %10.4e %s \n" % (j,dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])

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
            stri_p += "".join([ "%-10s %-10s %4s \n" % (res[pair[0][0][e][0]],res[pair[0][0][e][1]],pair[0][1][e]) for e in range(len(pair[0][0]))])
            
            stri_s += "# PDB %s \n" % args.pdbs[i].split("/")[-1]
            stri_s += "".join([ "%-10s %-10s %4s \n" % (res[st[0][0][e][0]],res[st[0][0][e][1]],st[0][1][e]) for e in range(len(st[0][0]))])

            if(args.dotbr):
                dotbr = bb.dot_bracket(pair,res)
                stri_dot += "%-20s %s\n" %(args.pdbs[i].split("/")[-1],dotbr[0])

            
    else:
        st,pair,res = bb.annotate(args.trj,topology=args.top)
        if(args.dotbr):
            dotbr = bb.dot_bracket(pair,res)
            stri_dot += "".join(["%-10d %s\n" %(k,dotbr[k]) for k in range(len(pair))])
            
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

##################### DUMP #######################

def dump(args):

    assert args.dumpR or args.dumpG, "# ERROR. choose --dumpR and/or --dumpR"
    

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
                angles_b,rr = bb.pucker_angles(args.pdbs[i],residues=args.res)
                stri_b += "".join(["%-12s %s \n" % (rr[e], "".join([" %11.3e" % angles_b[0,e,k] for k in range(angles_b.shape[2])])) for e in range(angles_b.shape[1])])
        else:
            
            angles_b,rr = bb.pucker_angles(args.trj,topology=args.top,residues=args.res)
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





####################### MAIN #########################

def main():

    def filename(args):
        # create output filename
        if(args.name == None):
            outfile = 'outfile.' + args.subparser_name 
        else:
            outfile = args.name + "." + args.subparser_name
        args.name = outfile
        print("# your output will be written to files with prefix %s" % outfile)

    # Parse options
    args = parse()

    # create filename
    if(args.subparser_name!="SNIPPET"):
        outfile = filename(args)

    # check
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
               'DUMP' : dump,\
               'TORSION':torsion,\
               'JCOUPLING':couplings,\
               'ENM':enm,\
               'SNIPPET' : snippet}

    options[args.subparser_name](args)
    
####################### MAIN ########################    


if __name__ == "__main__":
    main()
