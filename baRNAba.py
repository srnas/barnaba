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

import cProfile

import sys
import time

if sys.version_info[0] != 2 or sys.version_info[1] != 7:
    sys.stderr.write('# Python 2.7 is required. Aborting \n')
    sys.exit(1)
else:
    import argparse

def parse():

    parser = argparse.ArgumentParser(description='This is baRNAba')
    parser.add_argument("--name", dest="name",help="Job ID",default='',required=False)
    parser.add_argument("--res-mode", dest="res_mode",help="set mode: R:RNA,D:DNA,P:PROTEIN",default='R',required=False)

    subparsers = parser.add_subparsers(title="Subcommands",dest='subparser_name')

    # TESTED XTC
    parser_a = subparsers.add_parser('ERMSD', help='calculate ERMSD')
    parser_a.add_argument("--pdb", dest="reference",help="Reference PDB file",required=True)
    parser_a.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_a.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_a.add_argument("--per-res", dest="perres",help="Print per-residue eERMSD",action='store_true',default=False)
    parser_a.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)
    
    # TESTED XTC
    parser_b = subparsers.add_parser('ESCORE', help='Calculate Escore')
    parser_b.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=False)
    parser_b.add_argument("--force-field", dest="ff",help="PDB force field file",default='',required=True)
    parser_b.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=1.58,type=float)    
    parser_b.add_argument("--type", dest="type",default='standard',choices=['standard','beta'],
                              help='Type of ESCORE calculation (default=standard), beta not implemented')
    parser_b.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)

    # TESTED XTC
    parser_c = subparsers.add_parser('SS_MOTIF', help='Search single stranded (hairpin loop) RNA motif ')
    parser_c.add_argument("--query", dest="reference",help="Reference PDB file",required=True)
    parser_c.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_c.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_c.add_argument("--treshold", dest="treshold",help="ERMSD treshold",default=0.7,type=float)    
    parser_c.add_argument("--bulges", dest="bulges",help="Number of allowed bulged nucleotides",default=0,type=int)   
    parser_c.add_argument("--sequence", dest="seq",help="Sequence Accepts ACGU/NRY/ format. Default = any",required=False,default=None)
    parser_c.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)
    parser_c.add_argument("--dumpPDB", dest="dump_pdb",help="Write pdb files",action='store_true',default=False)

    parser_d = subparsers.add_parser('DS_MOTIF', help='Search RNA Double stranded (internal loop) motifs')
    parser_d.add_argument("--query", dest="reference",help="Reference PDB file",required=True)
    parser_d.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_d.add_argument("--l1", dest="l1",help="Length of first strand",required=True,type=int)    
    parser_d.add_argument("--l2", dest="l2",help="Lenght of second strand",required=True,type=int)    
    parser_d.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_d.add_argument("--treshold", dest="treshold",help="ERMSD treshold",default=0.7,type=float)    
    parser_d.add_argument("--bulges", dest="bulges",help="Number of allowed bulged nucleotides per strand)",default=0,type=int)    
    parser_d.add_argument("--sequence", dest="seq",help="Sequence Accepts ACGU/NRY/ format. Default = any",required=False,default=None)
    parser_d.add_argument("--dumpPDB", dest="dump_pdb",help="Write pdb files",action='store_true',default=False)
    parser_d.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)
    
    # TESTED XTC
    parser_e = subparsers.add_parser('ANNOTATE', help='Annotate RNA structure')
    parser_e.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_e.add_argument("--pymol", dest="pymol",help="Write script.pml file for coloring",action='store_true',default=False)
    parser_e.add_argument("--hread", dest="hread",help="make output human-readable",action='store_true',default=True)
    parser_e.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)
    
    # TESTED XTC
    parser_f = subparsers.add_parser('DUMP', help='DUMP structural parameters')
    parser_f.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_f.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_f.add_argument("--dumpG", dest="dumpG",help="Write G vectors on .gvec file",action='store_true',default=False)
    parser_f.add_argument("--dumpR", dest="dumpR",help="Write R vectors on .rvec file",action='store_true',default=False)
    parser_f.add_argument("--dumpP", dest="dumpP",help="Write P vectors on .pvec file",action='store_true',default=False)
    parser_f.add_argument("--hread", dest="read",help="make output human-readable",action='store_true',default=False)
    parser_f.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)

    # TESTED XTC
    parser_g = subparsers.add_parser('SNIPPET', help='SPLIT structure in multiple PDB')
    parser_g.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_g.add_argument("--seq", dest="seq",help="Sequence type. Accepts ACGU/NRY format",required=True)
    parser_g.add_argument("--cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)

    # TESTED XTC
    parser_h = subparsers.add_parser('TORSION', help='Calculate dihedral angles')
    parser_h.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_h.add_argument("--hread", dest="hread",help="make output human-readable",action='store_true',default=False)
    parser_h.add_argument("--backbone",dest="bb",help="calculate torsion backbone angles",action="store_true",default=False)
    parser_h.add_argument("--pucker",dest="pucker",help="calculate pucker angles",action="store_true",default=False)
    parser_h.add_argument("--jcouplings",dest="jcoupling",help="calculate torsion backbone angles",action="store_true",default=False)
    parser_h.add_argument("--raw",dest="raw",help="print raw angles for j3",action="store_true",default=False)
    parser_h.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)
                          
    # TESTED XTC
    parser_j = subparsers.add_parser('ENM', help='Calculate ENM')
    parser_j.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_j.add_argument("--cutoff", dest="cutoff",help="Cutoff distance in Angstrom (default=9)",default=9.0,type=float)
    parser_j.add_argument("--type", dest="type",default='SBP',choices=['P','S','B','SBP','AA'], help='Type of ENM (default=SBP)')    
    parser_j.add_argument("--ntop", dest="ntop",help="Number of top eigenvectors to write (default=10)",default=10,type=int)
    parser_j.add_argument("--zmodes", dest="zmodes",help="Write modes corresponding to zero eigenvalues",action='store_true',default=False)

    parser_k = subparsers.add_parser('NOE', help='Calculate NOE')
    parser_k.add_argument("-f", dest="files",help="PDB file(s)",nargs="+",default='',required=True)
    parser_k.add_argument("--xtc",dest="xtc",help="XTC gromacs trajectory",required=False,default=None)
    parser_k.add_argument("--cutoff", dest="cutoff",help="Print to file only average distance below cutoff",default=7.0,type=float)  
    
    args = parser.parse_args()

    return args


####################### ERMSD ########################

def ermsd(args):

    import ermsd    
    ermsd.ermsd(args)


####################### SCORING ########################
        
def score(args):

    import score
    score.score(args)


####################### MOTIF #########################

def ss_motif(args):

    import ss_motif
    ss_motif.ss_motif(args)


def ds_motif(args):

    import ds_motif
    ds_motif.ds_motif(args)
    
##################### ANNOTATE #######################

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

####################  NOE  #######################

def noe(args):
    import noe
    noe.noe(args)

####################### MAIN #########################

def main():

    def filename(args):
        # create output filename
        if(args.name == ''):
            w = time.strftime("%j.%m.%M.")
            outfile = w + args.subparser_name + ".rna"
            args.name = outfile
            print "# No name specified. Your output will be written to",outfile
        else:
            outfile = args.name + "." + args.subparser_name + ".rna"
            print "# Your output will be written to",outfile
            args.name = outfile

    
    # check at startup!
    try:
        import numpy
    except ImportError:
        sys.stderr.write('# Numpy is not installed \n')
        sys.exit(1)

    try:
        import scipy
    except ImportError:
        sys.stderr.write('# Scipy is not installed \n')
        sys.exit(1)

    # Parse options
    args = parse()

    # check pdb file
    for f in args.files:
        assert f[-4:]==".pdb", "# Error: PDB format only"


    # create filename
    outfile = filename(args)

    # call appropriate function
    options = {'ERMSD' : ermsd,'ESCORE' : score,'SS_MOTIF' : ss_motif,'DS_MOTIF' : ds_motif,\
                   'ANNOTATE' : annotate,'DUMP' : dump,'TORSION':torsion,'ENM':enm,\
                   'SNIPPET' : snippet,"NOE":noe}

    options[args.subparser_name](args)
    
####################### MAIN ########################    


if __name__ == "__main__":
    main()
