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

import numpy as np
import mdtraj as md
import btools as bt
import definitions as definitions
from scipy.spatial import distance

####################### MOTIF #########################


def ss_motif(args):

    print "# Finding 3D Single Strand Motifs..."

    # load reference structure
    ref_pdb = md.load(args.reference)
    ref_len = ref_pdb.n_residues

    # calculate gmat
    ref_mat = (bt.pdb2gmat(ref_pdb,args.cutoff)).reshape(-1)

    # get indeces for superposition 
    ref_sup_idx = bt.get_sup_idx(ref_pdb.topology)

    # assign query sequence (default N*ll)
    query="N"*ref_len
    if(args.seq!=None): query = args.seq
    assert len(query)==ref_len, "# FATAL: query structure and sequence length mismatch!"

    fh = open(args.name,'a')
    if(args.pdbs!=None):
        files = args.pdbs
        string = "#%s %8s %s %s \n"  % ("Index","Distance","Sequence","File")
    else:
        files = [args.trj]
        string = "#%s %8s %s %s \n"  % ("Index","Distance","Sequence","Time")
    fh.write(string)

    # loop over files
    for i in xrange(0,len(files)):

        name_pref = files[i][0:-4].split("/")[-1]
                
        if(args.pdbs!=None):
            cur_pdb = md.load_pdb(files[i])
        else:
            top=args.top
            cur_pdb = md.load_frame(files[i],0,top=top)

        if(cur_pdb.n_residues<ref_pdb.n_residues): continue

        cur_idx = bt.get_lcs_idx(cur_pdb.topology)
        ii = 0
        
        # get list of residues
        rna_residues = [cur_pdb.topology.atom(at).residue  for at in cur_idx[0]]
        # get sequence
        rna_seq = [definitions.residue_dict[rr.name]  for rr in rna_residues]
        rna_seq_id = ["%s_%d_%d" % (rr.name,rr.resSeq,rr.chain.index)  for rr in rna_residues]
        res_idxs_full = [rr.index  for rr in rna_residues]
        
       # find indeces in rna_seq matching query
        res_idxs = bt.get_idx(rna_seq,query,args.bulges)
        if(len(res_idxs)==0): return 0
        
        if(args.pdbs!=None):
            
            # coordinates for LCS calculations
            coords1 = cur_pdb.xyz[0,cur_idx[0]]
            coords2 = cur_pdb.xyz[0,cur_idx[1]]
            coords3 = cur_pdb.xyz[0,cur_idx[2]]
            # get gmats
            gmats = [[bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in res_idxs]]

            # calculate distances
            dists = [distance.cdist([ref_mat],el)/np.sqrt(ref_len) for el in gmats]
            below_t = [(dd<args.treshold).nonzero()[1] for dd in dists]
            
            for t in range(len(below_t)):
                for it in below_t[t]:
                    idxs = res_idxs[it]
                    seq = " ".join([rna_seq_id[p] for p in idxs ])
                    fh.write("%4d %8.5f %80s %50s \n" % (ii,dists[t][0,it],seq,files[i]))
                    ii += 1

                    if(args.dump==False): continue
                    
                    resi_full = [cur_pdb.topology.residue(res_idxs_full[y]) for y in idxs]
                    tmp_atoms0 = [[atom.index for atom in resi.atoms] for resi in resi_full]
                    # flatten
                    tmp_atoms = [val for sublist in tmp_atoms0 for val in sublist]
                    trj_slice = cur_pdb.atom_slice(tmp_atoms)
                    
                    # get atoms for alignment
                    cur_sup_idx = bt.get_sup_idx(trj_slice.topology)
                    
                    # align and save
                    name = name_pref + "_" + str(ii).zfill(4) + ".pdb"                    
                    if(len(cur_sup_idx)==len(ref_sup_idx)):
                        trj_slice.superpose(ref_pdb,atom_indices=cur_sup_idx,ref_atom_indices=ref_sup_idx)
                    trj_slice.save(name)

                        
        else:

            ii = 0
            # trajectory: split in chunks
            for chunk in md.iterload(files[i], chunk=100,top=top):
                # coordinates for LCS calculations
                coords1 = chunk.xyz[:,cur_idx[0]]
                coords2 = chunk.xyz[:,cur_idx[1]]
                coords3 = chunk.xyz[:,cur_idx[2]]
                
                gmats = [[bt.get_gmat(coords1[t,idx],coords2[t,idx],coords3[t,idx],args.cutoff).reshape(-1) for idx in res_idxs] for t in range(len(chunk))]

                dists = [distance.cdist([ref_mat],el)/np.sqrt(ref_len) for el in gmats]
                below_t = [(dd<args.treshold).nonzero()[1] for dd in dists]
                for t in range(len(below_t)):
                    for it in below_t[t]:
                        seq = " ".join([rna_seq_id[p] for p in res_idxs[it] ])
                        fh.write("%4d %8.5f %80s %10f\n" % (ii,dists[t][0,it],seq,chunk[t].time))
                        ii += 1

                        if(args.dump==False): continue
                        idxs = res_idxs[it]
                    
                        resi_full = [cur_pdb.topology.residue(res_idxs_full[y]) for y in idxs]
                        tmp_atoms0 = [[atom.index for atom in resi.atoms] for resi in resi_full]
                        # flatten
                        tmp_atoms = [val for sublist in tmp_atoms0 for val in sublist]
                        trj_slice = cur_pdb.atom_slice(tmp_atoms)
                        
                        # get atoms for alignment
                        cur_sup_idx = bt.get_sup_idx(trj_slice.topology)
                        
                        # align and save
                        name = name_pref + "_" + str(ii).zfill(4) + ".pdb"                    
                        if(len(cur_sup_idx)==len(ref_sup_idx)):
                            trj_slice.superpose(ref_pdb,atom_indices=cur_sup_idx,ref_atom_indices=ref_sup_idx)
                        trj_slice.save(name)
                    
    return 0




