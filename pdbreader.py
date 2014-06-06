import sys
import os
import numpy as N
import time

rnas = ['U','rU','RU','RU5','RU3','C','rC','RC','RC5','RC3',\
            'G','rG','RG','RG5','RG3','A','rA','RA','RA5','RA3']       
pyr = ['U','rU','RU','RU5','RU3','C','rC','RC','RC5','RC3']
pur = ['G','rG','RG','RG5','RG3','A','rA','RA','RA5','RA3']
types = ['C4','C2','C6']


def write_args(args,fh):
    fh.write("# This is a baRNAba run. Timestamp: " + time.strftime("%c\n"))
    for k in args.__dict__:
        if(str(k) != 'files'):
            s = "# " + str(k) + " " + str(args.__dict__[k]) + "\n"
            fh.write(s)


def get_coord(file):

    def rearrange(raw_coords,raw_sequence):
        # do some extra checks
        if(len(raw_coords)%3!=0):
            counts = [1]
            ss = []
            for i in range(1,len(raw_sequence)):
                if(raw_sequence[i][0]==raw_sequence[i-1][0]):
                    counts[-1] += 1
                else:
                    counts.append(1)
                    ss.append(raw_sequence[i][0])
            for i in range(len(counts)):
                if(counts[i] != 3):
                    print "# Only",counts[i],"atoms in residue",ss[i-1]
            print "# FATAL ERROR. Check your PDB file"
            sys.exit(1)

        ll = len(raw_coords)/3
        if(ll==0):
            print "# Fatal error. No RNA coordinates in PDB. Exiting"
            sys.exit(1)
        if(ll==1):
            print "# Warning: only 1 RNA residue in PDB."

            
        # rearrange array (swap C4/C6)
        coords = N.zeros((ll,3,3))
        sequence = []
        for i in range(ll):
            sequence.append(raw_sequence[3*i][0])
            for j in range(3):
                j1 = get_idx(raw_sequence[3*i+j])
                for k in range(3):
                    coords[i,j,k] = raw_coords[3*i+j1][k]
        return coords,sequence,len(sequence)


    def get_idx(v):

        atm = v[1]
        res = v[0].split("_")[1]
        if(atm == 'C2'):
            return 0

        if(atm == 'C4' and res in pyr):
            return 1
        if(atm == 'C4' and res in pur):
            return 2
        if(atm == 'C6' and res in pyr):
            return 2
        if(atm == 'C6' and res in pur):
            return 1
        return -3


    print "# Reading file ",file,
    raw_coords = []
    raw_sequence = []
    coords = []
    sequences = []
    stat = []

    fh = open(file)
    for line in fh:           
        s = line.split()
        if(len(s)==0):
            continue
        # check line
        if (s[0]=="ENDMDL"):
            c,s,n = rearrange(raw_coords,raw_sequence)
            stat.append(n)
            coords.append(c)
            sequences.append(s)
            raw_coords = []
            raw_sequence = []

        if(s[0] != "ATOM"):
            continue
        # check residue type
        REST = line[17:20].strip()
        if(REST not in rnas):
            continue
        if(len(REST)>1):
            REST = REST[1]

        # check atom type
        ATMT = line[12:16].strip()
        if(ATMT not in types):
            continue


        # check alternate location
        ALT = line[16:17].strip()
        if(ALT != "" and ALT != "A"):
            continue

        RESN = line[22:26].strip()
        CHAIN = line[21:22].strip()
        RES_ID = RESN + "_" + REST + "_" + CHAIN
        raw_sequence.append([RES_ID,ATMT])
        raw_coords.append([float(line[30:38]),float(line[38:46]),float(line[46:53])])

    if(len(coords)==0):
        c,s,n = rearrange(raw_coords,raw_sequence)
        coords.append(c)
        stat.append(n)
        sequences.append(s)

    print "...",len(coords),"model(s) read. "

    return coords,sequences



def pdb_string(atom_nr,atom_name,res_name,chain,res_nr,x,y,z,bfac,occ,tt):
    string = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f           %2s%2s \n" \
        % ("ATOM",atom_nr,atom_name," ",res_name,chain,res_nr," ",x,y,z,bfac,occ,tt,' ')
    return string



def create_ndx(pdb,ndx_file):
    
    array = []
    fh = open(pdb)
    for line in fh:
        # check line
        s = line.split()
        if (len(s) < 9):
            continue
        if(s[0] != "ATOM"):
            continue
        # check residue type
        REST = line[17:20].strip()
        if(REST not in rnas):
            continue
        if(len(REST)>1):
            REST = REST[1]
        # check atom type
        ATMT = line[12:16].strip()
        if(ATMT not in types):
            continue
        RESN = int(line[22:26].strip())
        ATMN = line[6:11].strip()
        array.append([RESN,REST,ATMT,ATMN])
    fh.close()

    ll = len(array)/3
    ii = -1*N.ones((int(ll),3))
    for el in array:
        j = el[0] - array[0][0]
        if(el[2] == 'C2'):
            ii[j][0] = el[3] 
        if(el[2] == 'C4' and el[1] in pyr):
            ii[j][1] = el[3] 
        if(el[2] == 'C4' and el[1] in pur):
            ii[j][2] = el[3] 
        if(el[2] == 'C6' and el[1] in pyr):
            ii[j][2] = el[3] 
        if(el[2] == 'C6' and el[1] in pur):
            ii[j][1] = el[3] 
    
    idx = [i for ss in ii for i in ss]
    string  = '[ C2_C4_C6 ] \n'
    for el in idx:
        if(el==-1):
            print "# Fatal Error: check the pdb file"
            sys.exit(1)
        string += str(int(el)) + " "
    string += "\n"
    fh_out = open(ndx_file,'w')
    fh_out.write(string)
    fh_out.close()

def xtc2pdb(xtc,args):


    # create index file
    ndx_file = 'tmp_index.ndx'
    create_ndx(args.pdb,ndx_file)
    
    # print concatenated pdb
    output = (args.name).split(".rna")[0] + "_samples.pdb"
    cmd = "trjconv" + args.gmx + \
        " -f " +  xtc + \
        " -o " + output + \
        " -n " + ndx_file + \
        " -s " + args.pdb + \
        " -skip " + str(args.skip) + " &> gmx.tmp"

    print "# You provided an XTC trajectory."
    print "# Make ABSOLUTELY sure that: "
    print "# 1) PBC were removed "
    print "# 2) the indexing in the XTC and in the reference PDB file coincide"
    print "# Converting trajectory: ", cmd
    os.system(cmd)
    print "# Coordinates of relevant atoms written to ", output
    return [output]

