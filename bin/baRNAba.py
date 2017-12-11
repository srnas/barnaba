#!/usr/bin/env python2.7

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


import sys
import glob
import os
import argparse
import barnaba as bb
    
def parse():

    parser = argparse.ArgumentParser(description='This is baRNAba')
    subparsers = parser.add_subparsers(title="Subcommands",dest='subparser_name')

    # ERMSD PARSER
    parser_01 = subparsers.add_parser('ERMSD', help='calculate eRMSD')
    parser_01.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_01.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",required=False,default=None)
    parser_01.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_01.add_argument("--top", dest="top",help="Topology file",required=False)
    
    parser_01.add_argument("--ref", dest="reference",help="Reference PDB file",required=True)
    parser_01.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)

    
    #############
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
    parser_05.add_argument("--hread", dest="hread",help="make output human-readable",action='store_true',default=False)

    
    # DUMP
    parser_06 = subparsers.add_parser('DUMP', help='DUMP structural parameters')
    parser_06.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_06.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_06.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_06.add_argument("--top", dest="top",help="Topology file",required=False)
    parser_06.add_argument("--stride", dest="stride",help="Stride",required=False,default=1,type=int)

    parser_06.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_06.add_argument("--dumpG", dest="dumpG",help="Write G vectors on .gvec file",action='store_true',default=False)
    parser_06.add_argument("--dumpR", dest="dumpR",help="Write R vectors on .rvec file",action='store_true',default=False)
    parser_06.add_argument("--dumpP", dest="dumpP",help="Write P vectors on .pvec file",action='store_true',default=False)
    parser_06.add_argument("--atom", dest="atom",help="Specify atom type (default = P)",default="P",required=False)
    parser_06.add_argument("--hread", dest="read",help="make output human-readable",action='store_true',default=False)
    
    # CREATE FRAGMENTS
    parser_07 = subparsers.add_parser('SNIPPET', help='SPLIT structure in multiple PDB')
    parser_07.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_07.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=True)
    parser_07.add_argument("--seq", dest="seq",help="Sequence type. Accepts ACGU/NRY format",required=True)


    # ###########################################
    # CALCULATE TORSION ANGLES
    parser_08 = subparsers.add_parser('TORSION', help='Calculate dihedral angles')
    parser_08.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_08.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_08.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_08.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_08.add_argument("--hread", dest="hread",help="make output human-readable",action='store_true',default=False)

    # CALCULATE TORSION ANGLES - J COUPLINGS OK
    parser_09 = subparsers.add_parser('JCOUPLING', help='Calculate J couplings.')
    parser_09.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_09.add_argument("--pdb", dest="pdbs",help="PDB file(s)",nargs="+",default=None,required=False)
    parser_09.add_argument("--trj", dest="trj",help="Trajectory",required=False,default=None)
    parser_09.add_argument("--top", dest="top",help="Topology file",required=False)

    parser_09.add_argument("--hread", dest="hread",help="make output human-readable",action='store_true',default=False)
    parser_09.add_argument("--raw",dest="raw",help="print raw angles for j3",action="store_true",default=False)
                          
    # CALCULATE ENM
    parser_10 = subparsers.add_parser('ENM', help='Calculate ENM')
    parser_10.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_10.add_argument("--pdb", dest="pdbs",help="PDB file",default=None,required=True)

    parser_10.add_argument("--cutoff", dest="cutoff",help="Cutoff distance in nm (default=0.9)",default=0.9,type=float)
    parser_10.add_argument("--type", dest="type",default='SBP',choices=['P','S','B','SBP','AA'], help='Type of ENM (default=SBP)')    
    parser_10.add_argument("--ntop", dest="ntop",help="Number of top eigenvectors to write (default=10)",default=10,type=int)
    parser_10.add_argument("--zmodes", dest="zmodes",help="Write modes corresponding to zero eigenvalues",action='store_true',default=False)
    parser_10.add_argument("--protein", dest="protein",help="Consider Proteins as well",action='store_true',default=False)


    # SMM
    parser_12 = subparsers.add_parser('SMM', help="Stop motion modeling")
    parser_12.add_argument("-o", dest="name",help="output_name",default=None,required=False)

    parser_12.add_argument("--gvec", dest="gvec",help="GVEC",required=True)
    parser_12.add_argument("--weight", dest="weight",help="optional weights for GVEC",required=False,default=None)
    parser_12.add_argument("--temp", dest="temp",help="Set temperature. This will interpret weights as free energies",required=False,default=0.0,type=float)
    parser_12.add_argument("--ref", dest="ref",help="GVEC of source and sink reference structure",required=True)

    parser_12.add_argument("--eps", dest="eps",help="eps",default=0.2,type=float)
    parser_12.add_argument("--epsc", dest="epsc",help="eps",default=0.7,type=float)
    parser_12.add_argument("--ncluster", dest="ncluster",required=False,default=20,type=int)
    parser_12.add_argument("--ntraj", dest="ntraj",required=False,default=0,type=int)
    parser_12.add_argument("--dist", dest="dist",help="distance for defining source and sink",default=0.3,type=float)
    parser_12.add_argument("--projection", dest="mode",required=False,default="kernel",choices=["kernel","PCA"])
    
    # CLUSTERING
    parser_13 = subparsers.add_parser('CLUSTER', help="Stop motion modeling")
    parser_13.add_argument("-o", dest="name",help="output_name",default=None,required=False)
    parser_13.add_argument("--gvec", dest="gvec",help="GVEC",required=True)
   
    parser_13.add_argument("--ncluster", dest="ncluster",required=False,default=20,type=int)
    #parser_13.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_13.add_argument("--proj", dest="mode",required=False,default="kernel",choices=["kernel","PCA"])
    parser_13.add_argument("--eps_p", dest="eps_p",help="eps for projection (used in kernel only)",default=0.8,type=float)
    parser_13.add_argument("--alg", dest="alg",required=False,default="kernel",choices=["kernel","DBSCAN"])
    parser_13.add_argument("--eps", dest="eps",help="epsilon used for clustering. default=0.7",default=0.7,type=float)
    parser_13.add_argument("--mins", dest="mins",help="minsample (for DBSCAN only), Default=50",default=50,type=float)
    
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
    
def rmsd(args):


    if(args.top==None):
        if(args.dump==True):
            dd = [bb.rmsd(args.reference,pdb,out="%s_%06d.pdb" % (out,i)) for i,pdb in enumerate(args.pdbs)]
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



def ss_motif(args):

    out = None
    if(args.dump==True):out = args.name

    stri = "# %s \n" % (" ".join(sys.argv[:]))
    if(args.top==None):
        stri += "#%20s %10s %s \n" % ("PDB","eRMSD","Sequence")
        for i in range(len(args.pdbs)):
            dd = bb.ss_motif(args.query,args.pdbs[i],out=out,bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
            stri += " ".join([" %20s %10.4e %s \n" % (args.pdbs[i].split("/")[-1],dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])
    else:
        dd = bb.ss_motif(args.query,args.trj,topology=args.top,out=out,bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
        stri += " ".join([" %20d %10.4e %s \n" % (j,dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])

    fh = open(args.name + ".out",'w')
    fh.write(stri)
    fh.close()


def ds_motif(args):

    out = None
    if(args.dump==True):out = args.name

    stri = "# %s \n" % (" ".join(sys.argv[:]))
    if(args.top==None):
        stri += "#%20s %10s %s \n" % ("PDB","eRMSD","Sequence")
        for i in range(len(args.pdbs)):
            dd = bb.ds_motif(args.query,args.pdbs[i],out=out,l1=args.l1,l2=args.l2,\
                             bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
            stri += " ".join([" %20s %10.4e %s \n" % (args.pdbs[i].split("/")[-1],dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])
    else:
        stri += "#%20s %10s %s \n" % ("frame","eRMSD","Sequence")
        dd = bb.ss_motif(args.query,args.trj,topology=args.top,out=out,l1=args.l1,l2=args.l2,\
                         bulges=args.bulges,threshold=args.threshold,sequence=args.seq,cutoff=args.cutoff)
        stri += " ".join([" %20d %10.4e %s \n" % (j,dd[j][1],"-".join(dd[j][2])) for j in range(len(dd))])

    fh = open(args.name + ".out",'w')
    fh.write(stri)
    fh.close()


def annotate(args):

    import annotate
    annotate.annotate(args)

##################### ANNOTATE #######################

def dump(args):

    import dump
    dump.dump(args)

 ##################### ANNOTATE #######################
    
def torsion(args):

    import torsions
    torsions.torsions(args)

def couplings(args):
        
    import couplings
    couplings.couplings(args)


####################  SPLIT #######################

def split(args):
    import split
    split.split(args)


####################  ENM  #######################

def enm(args):
    import enm
    enm.enm(args)

####################  SNIPPET  #######################

def snippet(args):
    import snippet
    snippet.snippet(args)


####################  SMM  #######################

def smm(args):
    import smm
    smm.smm(args)


def cluster(args):
    import cluster
    cluster.cluster(args)


####################### MAIN #########################

def main():

    def filename(args):
        # create output filename
        if(args.name == None):
            outfile = 'outfile.' + args.subparser_name 
            ll = glob.glob(outfile + "*")
            args.name = outfile
            if(len(ll)!=0):
                sys.stderr.write('# Creating backup file:  %s \n' % (outfile + ".backup." + str(len(ll))))
                os.system("cp " + outfile + " " + outfile + ".backup." + str(len(ll)))
            print "# No name specified. Your output will be written to",outfile
        else:
            outfile = args.name + "." + args.subparser_name
            print "# Your output will be written to",outfile
            args.name = outfile
    

    # Parse options
    args = parse()

    # create filename
    outfile = filename(args)

    # check
    #if(args.subparser_name!="SMM" and args.subparser_name!="CLUSTER" and args.subparser_name!="SNIPPET"):#
    #
    #    if(args.pdbs==None):
    #        assert args.trj != None, "# Specify either pdbs (--pdb) or a trajectory file"
    #        assert args.top != None, "# Please provide a topology file"
    #    else:
    #        if(args.subparser_name!="ENM"):
    #            assert args.trj == None, "# Specify either pdbs (--pdb) or a trajectory file, not both"
    #        
    
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
               'SNIPPET' : snippet,\
               'SMM':smm,\
               'CLUSTER':cluster}

    options[args.subparser_name](args)
    
####################### MAIN ########################    


if __name__ == "__main__":
    main()
