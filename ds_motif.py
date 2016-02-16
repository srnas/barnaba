#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014-2016 Sandro Bottaro (sbottaro@sissa.it)

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


    
def ds_motif(args):

    # sanity checks
    print "# Finding 3D Double Stranded Motifs..."
    ref_pdb = md.load(args.reference)
    ref_len = ref_pdb.n_residues

    ref_len1 = args.l1
    ref_len2 = args.l2 
    assert(ref_len == ref_len1+ref_len2)
    
    if(args.seq==None):
        query1="N"*ref_len1
        query2="N"*ref_len2
    else:
        query1=(args.seq).split("%")[0]
        query2=(args.seq).split("%")[1]
        assert len(query1)==ref_len1, "# FATAL: query structure and sequence length mismatch!"
        assert len(query2)==ref_len2, "# FATAL: query structure and sequence length mismatch!"


    # get matrix and pattern of full structure
    ref_idx = bt.get_lcs_idx(ref_pdb.topology)
    c01 = ref_pdb.xyz[0,ref_idx[0]]
    c02 = ref_pdb.xyz[0,ref_idx[1]]
    c03 = ref_pdb.xyz[0,ref_idx[2]]
    ref_mat = (bt.get_gmat(c01,c02,c03,args.cutoff)).reshape(-1)
    # get indeces for superposition 
    ref_sup_idx = bt.get_sup_idx(ref_pdb.topology)

    # get matrix and pattern of first strand
    c11 =  ref_pdb.xyz[0,ref_idx[0][0:ref_len1]]
    c12 =  ref_pdb.xyz[0,ref_idx[1][0:ref_len1]]
    c13 =  ref_pdb.xyz[0,ref_idx[2][0:ref_len1]]
    ref_mat1 = (bt.get_gmat(c11,c12,c13,args.cutoff)).reshape(-1)
    pattern1= bt.get_pattern(query1)

    # get matrix and pattern of second strand
    c21 =  ref_pdb.xyz[0,ref_idx[0][ref_len1:]]
    c22 =  ref_pdb.xyz[0,ref_idx[1][ref_len1:]]
    c23 =  ref_pdb.xyz[0,ref_idx[2][ref_len1:]]
    ref_mat2 = (bt.get_gmat(c21,c22,c23,args.cutoff)).reshape(-1)
    pattern2= bt.get_pattern(query2)
            

    # calculate center of mass distances
    # this will be used to prune the search!
    ref_com1 = np.sum((c11+c12+c13)/3.,axis=0)/c11.shape[0]
    ref_com2 = np.sum((c21+c22+c23)/3.,axis=0)/c21.shape[0]
    dcom =  np.sqrt(np.sum((ref_com1-ref_com2)**2))

    # open file and write header
    fh = open(args.name,'w')
    fh.write("#%s %8s %s %s %s \n"  % ("Index","Distance","Sequence","File","Time"))
    
    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]


    for i in xrange(0,len(files)):

        name_pref = files[i][0:-4].split("/")[-1]
                
        if(args.pdbs!=None):
            cur_pdb = md.load_pdb(files[i])
        else:
            top=args.top
            cur_pdb = md.load_frame(files[i],0,top=top)

        if(cur_pdb.n_residues<ref_pdb.n_residues): continue

        cur_idx = bt.get_lcs_idx(cur_pdb.topology)

        # index for result count
        kk = 0

        rna_residues = [cur_pdb.topology.atom(at).residue  for at in cur_idx[0]]
        rna_seq_id = ["%s_%d_%d" % (rr.name,rr.resSeq,rr.chain.index)  for rr in rna_residues]
        rna_seq = [definitions.residue_dict[rr.name]  for rr in rna_residues]
        res_idxs_full = [rr.index  for rr in rna_residues]


        all_idx1 = bt.get_idx(rna_seq,query1,args.bulges)
        if(len(all_idx1)==0): continue
        all_idx2 = bt.get_idx(rna_seq,query2,args.bulges)
        if(len(all_idx2)==0): continue

        # pdbs
        if(args.pdbs!=None):
            
            # coordinates for LCS calculations
            coords1 = cur_pdb.xyz[0,cur_idx[0]]
            coords2 = cur_pdb.xyz[0,cur_idx[1]]
            coords3 = cur_pdb.xyz[0,cur_idx[2]]
            
            # get gmats and calculate distance for strand 1
            gmats1 = [bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in all_idx1]
            dists1 = distance.cdist([ref_mat1],gmats1)/np.sqrt(ref_len1)
            below_t1 = (dists1<args.treshold).nonzero()[1]
            if(len(below_t1) ==0): continue
            

            # get gmats and calculate distance for strand 2
            gmats2 = [bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in all_idx2]
            dists2 = distance.cdist([ref_mat2],gmats2)/np.sqrt(ref_len2)
            below_t2 = (dists2<args.treshold).nonzero()[1]
            if(len(below_t2) ==0): continue
            
            all_com = (coords1+coords2+coords3)/3.
            # calculate center of mass
            idx1 = [all_idx1[k] for k in below_t1]
            com1 = [np.sum(all_com[k],axis=0)/ref_len1 for k in idx1]
            if(len(com1) == 0): continue
            
            idx2 = [all_idx2[k] for k in below_t2]
            com2 = [np.sum(all_com[k],axis=0)/ref_len2 for k in idx2]
            if(len(com2) == 0): continue                

            # calculate distance between centers of mass 
            dmat = distance.cdist(com1,com2)
            c_idx = (dmat<1.5*dcom).nonzero()
            if(len(c_idx) == 0): continue
                    
            # get comboi indeces
            idx_combo_tmp = [idx1[c_idx[0][ii]] + idx2[c_idx[1][ii]] for ii in range(len(c_idx[0]))]
            # remove overlaps
            idx_combo = [el for el in idx_combo_tmp if(len(set(el)) == len(el))]
                    
            gmatsf = [bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in idx_combo]
            distsf = distance.cdist([ref_mat],gmatsf)/np.sqrt(ref_len)
            below_tf = (distsf<args.treshold).nonzero()[1]
                    
            for it in below_tf:
                idxs = idx_combo[it]
                seq = " ".join([rna_seq_id[p] for p in idxs ])
                fh.write("%4d %8.5f %80s %50s \n" % (kk,distsf[0,it],seq,files[i]))
                kk += 1

                if(args.dump==False): continue
                    
                resi_full = [cur_pdb.topology.residue(res_idxs_full[y]) for y in idxs]
                tmp_atoms0 = [[atom.index for atom in resi.atoms] for resi in resi_full]
                # flatten
                tmp_atoms = [val for sublist in tmp_atoms0 for val in sublist]
                trj_slice = cur_pdb.atom_slice(tmp_atoms)
                    
                # get atoms for alignment
                cur_sup_idx = bt.get_sup_idx(trj_slice.topology)
                    
                # align and save
                name = name_pref + "_" + str(kk).zfill(4) + ".pdb"                    
                if(len(cur_sup_idx)==len(ref_sup_idx)):
                    trj_slice.superpose(ref_pdb,atom_indices=cur_sup_idx,ref_atom_indices=ref_sup_idx)
                trj_slice.save(name)

        # trajectory
        else:

            for chunk in md.iterload(files[i], chunk=100,top=top):
                for j in range(len(chunk)):
                    # coordinates for LCS calculations
                    coords1 = chunk.xyz[j,cur_idx[0]]
                    coords2 = chunk.xyz[j,cur_idx[1]]
                    coords3 = chunk.xyz[j,cur_idx[2]]
                    
                    # get gmats and calculate distance for strand 1
                    gmats1 = [bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in all_idx1]
                    dists1 = distance.cdist([ref_mat1],gmats1)/np.sqrt(ref_len1)
                    below_t1 = (dists1<args.treshold).nonzero()[1]
                    if(len(below_t1) ==0): continue
            

                    # get gmats and calculate distance for strand 2
                    gmats2 = [bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in all_idx2]
                    dists2 = distance.cdist([ref_mat2],gmats2)/np.sqrt(ref_len2)
                    below_t2 = (dists2<args.treshold).nonzero()[1]
                    if(len(below_t2) ==0): continue
            
                    all_com = (coords1+coords2+coords3)/3.
                    # calculate center of mass
                    idx1 = [all_idx1[k] for k in below_t1]
                    com1 = [np.sum(all_com[k],axis=0)/ref_len1 for k in idx1]
                    if(len(com1) == 0): continue
            
                    idx2 = [all_idx2[k] for k in below_t2]
                    com2 = [np.sum(all_com[k],axis=0)/ref_len2 for k in idx2]
                    if(len(com2) == 0): continue                

                    # calculate distance between centers of mass 
                    dmat = distance.cdist(com1,com2)
                    c_idx = (dmat<1.5*dcom).nonzero()
                    if(len(c_idx) == 0): continue
                    
                    # get comboi indeces
                    idx_combo_tmp = [idx1[c_idx[0][ii]] + idx2[c_idx[1][ii]] for ii in range(len(c_idx[0]))]
                    # remove overlaps
                    idx_combo = [el for el in idx_combo_tmp if(len(set(el)) == len(el))]
                    
                    gmatsf = [bt.get_gmat(coords1[idx],coords2[idx],coords3[idx],args.cutoff).reshape(-1) for idx in idx_combo]
                    distsf = distance.cdist([ref_mat],gmatsf)/np.sqrt(ref_len)
                    below_tf = (distsf<args.treshold).nonzero()[1]
                    
                    for it in below_tf:
                        idxs = idx_combo[it]
                        seq = " ".join([rna_seq_id[p] for p in idxs ])
                        fh.write("%4d %8.5f %80s %10f \n" % (kk,distsf[0,it],seq,chunk[j].time))
                        kk += 1

                        if(args.dump==False): continue
                        
                        resi_full = [cur_pdb.topology.residue(res_idxs_full[y]) for y in idxs]
                        tmp_atoms0 = [[atom.index for atom in resi.atoms] for resi in resi_full]
                        # flatten
                        tmp_atoms = [val for sublist in tmp_atoms0 for val in sublist]
                        trj_slice = chunk[j].atom_slice(tmp_atoms)
                    
                        # get atoms for alignment
                        cur_sup_idx = bt.get_sup_idx(trj_slice.topology)
                    
                        # align and save
                        name = name_pref + "_" + str(kk).zfill(4) + ".pdb"                    
                        if(len(cur_sup_idx)==len(ref_sup_idx)):
                            trj_slice.superpose(ref_pdb,atom_indices=cur_sup_idx,ref_atom_indices=ref_sup_idx)
                        trj_slice.save(name)


    fh.close()


