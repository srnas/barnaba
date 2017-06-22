import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist
import sys
import nucleic
import definitions


def dsmotif_traj(ref,traj,l1,l2,treshold=0.9,cutoff=2.4,sequence=None,bulges=0,write=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = ref.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    ll = len(nn_ref.ok_residues)
    assert(ll==l1+l2)
    if(sequence==None):
        sequence1 = "N"*l1
        sequence2 = "N"*l2
    else:
        sequence1=(sequence).split("%")[0]
        sequence2=(sequence).split("%")[1]
        assert len(sequence1)==l1, "# FATAL: query structure and sequence length mismatch!"
        assert len(sequence2)==l2, "# FATAL: query structure and sequence length mismatch!"

        
    coords_ref = ref.xyz[0,nn_ref.indeces_lcs]
    ref_mat = nn_ref.get_gmat(coords_ref,cutoff).reshape(-1)

    coords_ref1 = ref.xyz[0,nn_ref.indeces_lcs[:,0:l1]]

    ref_mat1 = nn_ref.get_gmat(coords_ref1,cutoff).reshape(-1)
    coords_ref2 = ref.xyz[0,nn_ref.indeces_lcs[:,l1:]]
    ref_mat2 = nn_ref.get_gmat(coords_ref2,cutoff).reshape(-1)

   
    # calculate center of mass distances
    # this will be used to prune the search!
    ref_com1 = np.average(np.average(coords_ref1,axis=0),axis=0)
    ref_com2 = np.average(np.average(coords_ref2,axis=0),axis=0)
    dcom =  np.sqrt(np.sum((ref_com1-ref_com2)**2))
    
    # find indeces of residues according to sequence
    rna_seq = nn_traj.rna_seq_id
    all_idx1 = definitions.get_idx(rna_seq,sequence1,bulges)
    if(len(all_idx1)==0): return []
    all_idx2 = definitions.get_idx(rna_seq,sequence2,bulges)
    if(len(all_idx2)==0): return []

    # skip overlapping
    print all_idx1
    print all_idx2
    exit()
    #res_idxs = definitions.get_idx(rna_seq,sequence,bulges)
    #resname_idxs = [[nn_traj.rna_seq[l] for l in rr]  for rr  in res_idxs]
    #if(len(res_idxs)==0):
    #    return []

    lcs_idx = nn_traj.indeces_lcs
    
    results = []
    count = 1
    for i in xrange(traj.n_frames):

        # calculate center of mass distance
        
        gmats = []
        for j in res_idxs:
            coords_lcs = traj.xyz[i,lcs_idx[:,j]]
            gmats.append(nn_ref.get_gmat(coords_lcs,cutoff).reshape(-1))
        dd = cdist([ref_mat],gmats)/np.sqrt(ll)
        low = np.where(dd[0]<treshold)
        for k in low[0]:
            results.append([i,dd[0,k],resname_idxs[k]])

            # Write aligned PDB 
            if(write != None):
                pdb_out = "%s_%05d_%s_%d.pdb" % (write,count,resname_idxs[k][0],i)
                # slice trajectory
                tmp_atoms = []
                tmp_res =[]
                for r1 in res_idxs[k]:
                    tmp_atoms.extend([at.index for at in nn_traj.ok_residues[r1].atoms])
                    tmp_res.append([at for at in nn_traj.ok_residues[r1].atoms])
                traj_slice = traj[i].atom_slice(tmp_atoms)
                
                # align whatever is in common in the backbone
                idx_target = []
                idx_ref = []
                for res1,res2 in zip(nn_ref.ok_residues,traj_slice.topology.residues):
                    name2 = [at.name for at in res2.atoms if  at.name in definitions.bb_atoms]
                    for at in res1.atoms:
                        if at.name in definitions.bb_atoms:
                            if(at.name in name2):
                                idx_ref.append(at.index)
                                idx_target.append(((res2.atom(at.name)).index))
                                #idx_target.append(res2[name2.index(at.name)].index)
                traj_slice.superpose(ref,atom_indices=idx_target, ref_atom_indices=idx_ref)
                traj_slice.save(pdb_out)

            count += 1
    return results


def dsmotif(reference,target,l1,l2,treshold=0.9,cutoff=2.4,topology=None,sequence=None,bulges=0,write=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return dsmotif_traj(ref,traj,l1,l2,treshold=treshold,sequence=sequence,bulges=bulges,write=write)

