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
import mdtraj as md
import btools as bt
import numpy as np
import kde as kde

####################### SCORING ########################



def do_pdb(kernel,pdbs,cutoff):
    
    string = "#%29s %s \n " % ("Filename","ESCORE")
    
    for pdb in pdbs:
        cur_pdb = md.load_pdb(pdb)
        cur_idx = bt.get_lcs_idx(cur_pdb.topology)
        c1 = cur_pdb.xyz[0,cur_idx[0]]
        c2 = cur_pdb.xyz[0,cur_idx[1]]
        c3 = cur_pdb.xyz[0,cur_idx[2]]
        cur_mat = bt.get_mat_score(c1,c2,c3,cutoff+0.2)
        val = np.sum(kernel(10.0*cur_mat))
        string += '%30s %8.5f \n' % (pdb,val)
    return string

def do_traj(kernel,top,trj,cutoff):
    
    string = "#%9s %s \n " % ("time (ps)","ESCORE")

    # load first frame
    cur_pdb = md.load_frame(trj,0,top=top)
    cur_idx = bt.get_lcs_idx(cur_pdb.topology)
    
    # analyze trajectory in chunks of 100
    for chunk in md.iterload(trj, chunk=100,top=top):
        for j in range(len(chunk)):
            c1 = chunk.xyz[j,cur_idx[0]]
            c2 = chunk.xyz[j,cur_idx[1]]
            c3 = chunk.xyz[j,cur_idx[2]]
            cur_mat = bt.get_mat_score(c1,c2,c3,cutoff+0.2)
            val = np.sum(kernel(10.0*cur_mat))
            string += "%10.3f %8.5f \n" % (chunk.time[j],val)                        
    return string
    
def score(args):

    # calculate interaction matrix of the reference structure
    ref_pdb = md.load_pdb(args.reference)    
    ref_mat = bt.pdb2gmatscore(ref_pdb,args.cutoff)

    # set kernel density
    kernel = kde.gaussian_kde(10.0*ref_mat)
    
    # set kernel to 0.25 Angstrom
    kernel.set_bandwidth(0.25)

    print "# KDE computed. Bandwidth=",kernel.factor
    print "# Calculating ESCORE..."
    fh = open(args.name,'a')


    if(args.pdbs!=None):
        fh.write(do_pdb(kernel,args.pdbs,args.cutoff))
    else:
        fh.write(do_traj(kernel,args.top,args.trj,args.cutoff))
    fh.close()
    
    return 0

####################### SCORING ########################
