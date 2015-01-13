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
import time

if sys.version_info[0] != 2 or sys.version_info[1] != 7:
    print "# Python 2.7 is required. Aborting"
    sys.exit(1)
else:
    import argparse

def parse():

    parser = argparse.ArgumentParser(description='This is baRNAba')
    parser.add_argument("--name", dest="name",help="Job ID",default='',required=False)
    parser.add_argument("--trjconv", dest="gmx",help="Gromacs trjconv name. Necessary only when analyzing gro/xtc trajectories",default='',required=False)
    parser.add_argument("--skip", dest="skip",help="Skip frames in gro/xtc file",default='1',required=False)

 
    subparsers = parser.add_subparsers(title="Subcommands",dest='subparser_name')

    parser_a = subparsers.add_parser('ERMSD', help='calculate ERMSD')
    parser_a.add_argument("--pdb", dest="pdb",help="Reference PDB file",required=True)
    parser_a.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_a.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_a.add_argument("--type", dest="type",default='vector',choices=['scalar','vector'],
                              help='Type of ERMSD calculation')    
    parser_a.add_argument("--ermsf", dest="ermsf",help="Print per-residue ERMSD ",action='store_true',default=False)


    parser_b = subparsers.add_parser('ESCORE', help='Calculate Escore')
    parser_b.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=False)
    parser_b.add_argument("--force-field", dest="ff",help="PDB force field file",default='',required=True)
    parser_b.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=1.58,type=float)    
    parser_b.add_argument("--type", dest="type",default='standard',choices=['standard','beta'],
                              help='Type of ESCORE calculation (default=standard), beta not implemented')    
    parser_b.add_argument("--per-residue", dest="perres",help="Calculate per-residue score",action='store_true',default=False)

    parser_c = subparsers.add_parser('SS_MOTIF', help='Search single stranded (hairpin loop) RNA motif ')
    parser_c.add_argument("--query", dest="pdb",help="Reference PDB file",required=True)
    parser_c.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_c.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_c.add_argument("--treshold", dest="treshold",help="ERMSD treshold",default=0.8,type=float)    
    parser_c.add_argument("--bulges", dest="bulges",help="Number of allowed bulged nucleotides per strand)",default=0,type=int)   
    parser_c.add_argument("--type", dest="type",default='vector',choices=['scalar','vector'],
                              help='Type of ERMSD calculation (default=vector)')    

    parser_d = subparsers.add_parser('DS_MOTIF', help='Search RNA Double stranded (internal loop) motifs')
    parser_d.add_argument("--query", dest="pdb",help="Reference PDB file",required=True)
    parser_d.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_d.add_argument("--l1", dest="l1",help="Length of first strand",required=True,type=int)    
    parser_d.add_argument("--l2", dest="l2",help="Lenght of second strand",required=True,type=int)    
    parser_d.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_d.add_argument("--treshold", dest="treshold",help="ERMSD treshold",default=0.8,type=float)    
    parser_d.add_argument("--bulges", dest="bulges",help="Number of allowed bulged nucleotides per strand)",default=0,type=int)    
    parser_d.add_argument("--type", dest="type",default='vector',choices=['scalar','vector'],
                              help='Type of ERMSD calculation (default=scalar)')    

    parser_e = subparsers.add_parser('ANNOTATE', help='Annotate RNA structure')
    #parser_e.add_argument("--pdb", dest="pdb",help="Reference PDB file",required=False)
    parser_e.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=False)
    parser_e.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=1.58,type=float)    
    parser_e.add_argument("--compact", dest="compact",help="use compact format",action='store_true')    

    parser_f = subparsers.add_parser('DUMP', help='DUMP structural parameters')
    parser_f.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_f.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_f.add_argument("--type", dest="type",default='vector',choices=['scalar','vector'],
                              help='Type of calculation')    
    parser_f.add_argument("--dumpG", dest="dumpG",help="Write G vectors on .gvec file",action='store_true',default=False)
    parser_f.add_argument("--dumpR", dest="dumpR",help="Write R vectors on .rvec file",action='store_true',default=False)
    parser_f.add_argument("--dumpLCS", dest="dumpL",help="Write LCS on .lcs file",action='store_true',default=False)
    parser_f.add_argument("--hread", dest="read",help="make output human-readable",action='store_true',default=False)

    args = parser.parse_args()

    return args


####################### ERMSD ########################

def ermsd(args,files):

    import ermsd    
    ermsd.ermsd(args,files)


####################### SCORING ########################
        
def score(args,files):

    import score
    score.score(args,files)


####################### MOTIF #########################

def ss_motif(args,files):

    import ss_motif
    ss_motif.ss_motif(args,files)


def ds_motif(args,files):

    import ds_motif
    ds_motif.ds_motif(args,files)
    
##################### ANNOTATE #######################

def annotate(args,files):

    import annotate
    annotate.annotate(args,files)

##################### ANNOTATE #######################

def dump(args,files):

    import dump
    dump.dump(args,files)



####################### MAIN #########################

def main():

    def filename(args):
        # create output filename
        if(args.name == ''):
            #w = (args.pdb).split(".pdb")[0].split("/")[-1] + "."
            w = time.strftime("%j.%m.%M.")
            outfile = w + args.subparser_name + ".rna"
            args.name = outfile
            print "# No name specified. Your output will be written to",outfile
        else:
            outfile = args.name + "." + args.subparser_name + ".rna"
            print "# Your output will be written to",outfile
            args.name = outfile

    # check arguments and process file list
    def check(args):
        
        noref = ["ANNOTATE","ESCORE","DUMP"]
        if(args.subparser_name in noref):
            files = []
        else:
            files = [args.pdb]
        for f in args.files:
            extension = f[-4:]
            if(extension == '.pdb'):
                files.append(f)
            else:
                if(extension == '.xtc' or extension=='.gro'):
                    if(len(args.files) != 1):
                        print "# Fatal error. You cannot provide multiple XTC trajectories"
                        sys.exit(1)
                    else:
                        import pdbreader as pb
                        trj = pb.xtc2pdb(f,args)        
                        files.append(trj)
                else:
                    print "# Fatal error. Extension ",extension,"not recognized. Please check your -f file(s)"
                    sys.exit(1)
        return files

    # check at startup!
    try:
        import numpy
    except ImportError:
        print "# Numpy is not installed"
        sys.exit(1)

    try:
        import scipy
    except ImportError:
        print "# Scipy is not installed"
        sys.exit(1)

    # Parse options
    args = parse()

    # check xtc or pdb and return list
    files = check(args)

    # create filename
    outfile = filename(args)

    # call appropriate function
    options = {'ERMSD' : ermsd,\
                   'ESCORE' : score,\
                   'SS_MOTIF' : ss_motif,\
                   'DS_MOTIF' : ds_motif,\
                   'ANNOTATE' : annotate,\
                   'DUMP' : dump }
    options[args.subparser_name](args,files)
    
####################### MAIN ########################    


if __name__ == "__main__":
    main()
