import sys
import argparse

from scipy.spatial import distance
import numpy as N

import pdbreader as pb
import kde as k
import tools as t

import time

def parse():

    parser = argparse.ArgumentParser(description='This is baRNAba')
    parser.add_argument("-pdb", dest="pdb",help="Reference PDB file",required=True)
    parser.add_argument("-name", dest="name",help="Job ID",default='',required=False)
    parser.add_argument("-gmx", dest="gmx",help="GMX suffix",default='',required=False)
    parser.add_argument("-skip", dest="skip",help="Skip frames in XTC file",default='1',required=False)

 
    subparsers = parser.add_subparsers(title="Subcommands",dest='subparser_name')

    parser_a = subparsers.add_parser('ERMSD', help='calculate ERMSD')
    parser_a.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_a.add_argument("-cutoff", dest="cutoff",help="Ellipsoidal cutoff (default=2.4)",default=2.4,type=float)
    parser_a.add_argument("-type", dest="type",default='modulus',choices=['modulus','vector'],
                              help='Type of ERMSD calculation (default=modulus)')    
    parser_a.add_argument("-ermsf", dest="ermsf",help="Print per-residue ERMSD (to be implemented)",action='store_true')


    parser_b = subparsers.add_parser('ESCORE', help='Calculate Escore')
    parser_b.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=False)
    parser_b.add_argument("-force-field", dest="ff",help="PDB force field file",default='',required=True)
    parser_b.add_argument("-cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=1.58,type=float)    
    parser_b.add_argument("-type", dest="type",default='standard',choices=['standard','beta'],
                              help='Type of ESCORE calculation (default=standard), beta not implemented')    

    parser_c = subparsers.add_parser('SS_MOTIF', help='Search RNA motif ')
    parser_c.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_c.add_argument("-cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.0,type=float)    
    parser_c.add_argument("-treshold", dest="treshold",help="ERMSD treshold",default=0.35,type=float)    
    parser_c.add_argument("-type", dest="type",default='modulus',choices=['modulus','vector'],
                              help='Type of ERMSD calculation (default=modulus)')    

    parser_d = subparsers.add_parser('DS_MOTIF', help='Search RNA motif')
    parser_d.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=True)
    parser_d.add_argument("-l1", dest="l1",help="Length of first strand",default=-1,type=int)    
    parser_d.add_argument("-l2", dest="l2",help="Lenght of second strand",default=-1,type=int)    
    parser_d.add_argument("-cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=2.4,type=float)    
    parser_d.add_argument("-treshold", dest="treshold",help="ERMSD treshold",default=0.35,type=float)    
    parser_d.add_argument("-type", dest="type",default='modulus',choices=['modulus','vector'],
                              help='Type of ERMSD calculation (default=modulus)')    

    parser_e = subparsers.add_parser('ANNOTATE', help='Annotate RNA structure')
    parser_e.add_argument("-f", dest="files",help="PDB/XTC file(s)",nargs="+",default='',required=False)
    parser_e.add_argument("-cutoff", dest="cutoff",help="Ellipsoidal cutoff",default=1.58,type=float)    

    args = parser.parse_args()

    return args


def write_args(args,fh):
    fh.write("# This is a baRNAba run. Timestamp: " + time.strftime("%c\n"))
    for k in args.__dict__:
        if(str(k) != 'files'):
            s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
            fh.write(s)


######################## ERMSD #########################

def ermsd(args,files):
    

    fh = open(args.name,'w')

    print "# Calculating ERMSD..."
    write_args(args,fh)
    if(args.type=='modulus'):

        for ii,f in enumerate(files):
            atoms,sequence = pb.get_coord(f)

            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model. Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
                    mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                    ermsd = t.calc_dist_1d(ref_mat,mat)
                    string = '%8.5f %s.%i \n' % (ermsd,f,jj)
                    fh.write(string)
                

    if(args.type=='vector'):

        for ii,f in enumerate(files):
            atoms,sequence = pb.get_coord(f)

            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model. Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
                    mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                    ermsd = t.calc_dist_nd(ref_mat,mat)
                    string = '%8.5f %s.%i \n' % (ermsd,f,jj)
                    fh.write(string)

    fh.close()
    return 0

######################## ERMSD #########################


####################### SCORING ########################
        
def score(args,files):

    fh = open(args.name,'w')

    write_args(args,fh)

    ref_atoms,ref_sequence = pb.get_coord(args.ff)
    ref_lcs,ref_origo = t.coord2lcs(ref_atoms[0])
    ref_mat = t.lcs2mat_score(ref_lcs,ref_origo,args.cutoff)
    kernel = k.gaussian_kde(ref_mat)
    kernel.set_bandwidth(0.25)

    print "# KDE computed. Bandwidth=",kernel.factor
    print "# Calculating ESCORE..."

    for f in files:
        atoms,sequence = pb.get_coord(f)
        for ii,model in enumerate(atoms):
            lcs,origo = t.coord2lcs(model)

            # Cutoff is slightly augmented 
            mat = t.lcs2mat_score(lcs,origo,args.cutoff+0.2)
            val = kernel(mat)
            string = '%8.5f %s.%i \n' % (sum(val),f,ii)
            fh.write(string)

####################### SCORING ########################


####################### MOTIF #########################

def ss_motif(args,files):

    fh = open(args.name,'w')

    print "# Finding 3D motifs..."
    write_args(args,fh)

    
    if(args.type=='modulus'):

        for ii,f in enumerate(files):

            atoms,sequence = pb.get_coord(f)
            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model."
                    print "# Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                ll = ref_mat.shape[0]
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
      
                    mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                    for j in range(0,mat.shape[0]-ll+1):
                        red_mat = mat[j:j+ll,j:j+ll]
                        ermsd = t.calc_dist_1d(ref_mat,red_mat)
                        if(ermsd < args.treshold):
                            string = '%8.5f %s %i %s \n' % (ermsd,f,jj,sequence[jj][j:j+ll])
                            fh.write(string)


    if(args.type=='vector'):

        for ii,f in enumerate(files):

            atoms,sequence = pb.get_coord(f)
            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model."
                    print "# Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                ll = ref_mat.shape[0]
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
      
                    mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                    for j in range(0,mat.shape[0]-ll+1):
                        red_mat = mat[j:j+ll,j:j+ll]
                        ermsd = t.calc_dist_nd(ref_mat,red_mat)
                        if(ermsd < args.treshold):
                            string = '%8.5f %s %i %s \n' % (ermsd,f,jj,sequence[jj][j:j+ll])
                            fh.write(string)



def ds_motif(args,files):

    fh = open(args.name,'w')

    print "# Finding 3D motifs..."
    write_args(args,fh)
    l1 = args.l1
    l2 = args.l2

    if(args.type=='modulus'):

        for ii,f in enumerate(files):
            atoms,sequence = pb.get_coord(f)

            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model."
                    print "# Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                assert(ref_mat.shape[0] == args.l1+args.l2)
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
                    mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                    for j in range(0,mat.shape[0]-l1+1):
                        ends =  N.union1d(mat[j,:].nonzero()[0],mat[:,j].nonzero()[0])
                        #ends =  N.intersection1d(mat[j,:].nonzero()[0],mat[:,j].nonzero()[0])
                        for k in ends:
                            # standard case
                            if(k<l2-1):
                                continue
                            start1 = j
                            end1 = j+l1-1
                            start2 = k-l2+1
                            end2 = k
                            #print start1,end1,' ',start2,end2

                            # exclude if flanking pairs are not neighbors
                            if (mat[end1,start2] == 0 and mat[start2,end1] == 0):
                                continue
                            idx =  N.arange(start1,end1+1).tolist()
                            idx.extend(N.arange(start2,end2+1).tolist())
                            red_mat = mat[idx,:][:,idx]
                            #if(k<j):
                            #    red_mat = red_mat.T
                            ermsd = t.calc_dist_1d(ref_mat,red_mat)
                            #for aa in range(ref_mat.shape[0]):
                            #    for bb in range(ref_mat.shape[0]):
                            #        print aa+1,bb+1,ref_mat[aa,bb],red_mat[aa,bb],N.abs(ref_mat[aa,bb]-red_mat[aa,bb])

                            if(ermsd < args.treshold):
                                seq = ' - '
                                for el in sequence[jj][idx[0]:idx[l1-1]+1]:
                                    seq += (el + ' ')
                                seq += " - "    
                                for el in sequence[jj][idx[l1]:idx[-1]+1]:
                                    seq += (el + ' ')
                                    
                                string = '%8.5f %s %i %s \n' % (ermsd,f,jj,seq)
                                fh.write(string)
                                #print "#",string,
                                    
                                    #string = '%8.5f %s.%i' % (ermsd,f,jj)
                                    #fh.write(string)
                
####################### MOTIF ########################

##################### ANNOTATE #######################

def annotate(args,files):

    print "# Annotating RNA structures..."

    fh = open(args.name,'w')
    write_args(args,fh)

    ref_atoms,ref_sequence = pb.get_coord(files[0])
    idx = N.triu_indices(len(ref_sequence[0]),1)

    data = []    
    for f in files:
        atoms,sequence = pb.get_coord(f)
        for model,seq in zip(atoms,sequence):
            lcs,origo = t.coord2lcs(model)
            mat = t.lcs2mat_3d(lcs,origo,args.cutoff)
            int_mat = t.analyze_mat(mat,seq)
            data.append(int_mat[idx])


    data = N.array(data)

    # remove columns with no interactions and print to file
    idx1 = []
    header = "# "
    for k in range(data.shape[1]):
        if(any(data[:,k]) == 0):
            continue
        else:
            header += ref_sequence[0][idx[0][k]] + "/" + ref_sequence[0][idx[1][k]] + ":" + str(len(idx1)+2) + " "
            idx1.append(k)
    header += "\n"
    fh.write(header)

    # write data to file
    for k in range(data.shape[0]):
        string = "%6i" % k
        for l in idx1:
            string += "%4s " % (t.interactions[int(data[k,l])])
        fh.write(string + "\n")

    fh.close()
    return 0
            
##################### ANNOTATE #######################

####################### MAIN #########################

def main():

    def filename(args):
        # create output filename
        if(args.name == ''):
            w = (args.pdb).split(".pdb")[0].split("/")[-1] + "."
            w += time.strftime("%j.%M.")
            outfile = w + args.subparser_name + ".rna"
            args.name = outfile
            print "# No name specified. Your output will be written to",outfile
        else:
            outfile = args.name + "." + args.subparser_name + ".rna"
            print "# Your output will be written to",outfile
            args.name = outfile

    # check arguments and process file list
    def check(args):
        
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
                        files = pb.xtc2pdb(f,args)        
                else:
                    print "# Fatal error. Extension ",extension,"not recognized. Please check your -f file(s)"
                    sys.exit(1)
        return files

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
                   'ANNOTATE' : annotate}
    options[args.subparser_name](args,files)
    
####################### MAIN ########################    


if __name__ == "__main__":
    main()
