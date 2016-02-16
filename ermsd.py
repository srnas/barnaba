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

def calc_diff(ref,cur,pp):
    diff = (ref-cur)**2
    ll = ref.shape[0]
    val = np.sqrt(np.sum(diff)/ll)
    string  = "%10.6f " % (val)
    # calculate per-residue 
    if(pp):
        string += ";"
        per_res = np.sqrt(np.sum(np.sum(diff,axis=2),axis=1)/ll)
        for k in xrange(len(per_res)):
            string += " %10.6f " % (per_res[k])
    return string

def do_pdb(ref_mat,pdbs,cutoff,perres=False):
    
    string = "#%29s %s \n " % ("Filename","ERMSD")
    
    for pdb in pdbs:
        cur_pdb = md.load_pdb(pdb)
        assert(cur_pdb.n_residues==ref_mat.shape[0])
        # get indeces
        cur_idx = bt.get_lcs_idx(cur_pdb.topology)
        
        c1 = cur_pdb.xyz[0,cur_idx[0]]
        c2 = cur_pdb.xyz[0,cur_idx[1]]
        c3 = cur_pdb.xyz[0,cur_idx[2]]
        cur_mat = bt.get_gmat(c1,c2,c3,cutoff)
        ss = calc_diff(cur_mat,ref_mat,perres)
        string += "%30s %s \n" % (pdb,ss)
        
    return string


        
def do_traj(ref_mat,top,trj,cutoff,perres=False):

    string = "#%9s %s \n " % ("time (ps)","ERMSD")
    cur_pdb = md.load_frame(trj,0,top=top)
    assert(cur_pdb.n_residues==ref_mat.shape[0])

    # get indeces
    cur_idx = bt.get_lcs_idx(cur_pdb.topology)

    for chunk in md.iterload(trj, chunk=100,top=top):
        for j in range(len(chunk)):
            c1 = chunk.xyz[j,cur_idx[0]]
            c2 = chunk.xyz[j,cur_idx[1]]
            c3 = chunk.xyz[j,cur_idx[2]]
            cur_mat = bt.get_gmat(c1,c2,c3,cutoff)
            ss = calc_diff(cur_mat,ref_mat,perres)
            string += "%10.3f %s \n" % (chunk.time[j],ss)
    return string

def ermsd(args):
    
    print "# Calculating ERMSD..."

    # calculate interaction matrix of the reference structure
    
    ref_pdb = md.load(args.reference)
    ref_mat = bt.pdb2gmat(ref_pdb,args.cutoff)

    fh = open(args.name,'a')

    if(args.trj==None):
        fh.write(do_pdb(ref_mat,args.pdbs,args.cutoff,args.perres))
    else:
        fh.write(do_traj(ref_mat,args.top,args.trj,args.cutoff,args.perres))
    
    fh.close()
    return 0


