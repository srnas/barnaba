import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist
import tools
import sys
import nucleic
import definitions
import itertools 

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
    ref_mat = tools.calc_gmat(coords_ref,cutoff).reshape(-1)

    coords_ref1 = ref.xyz[0,nn_ref.indeces_lcs[:,0:l1]]

    ref_mat1 = tools.calc_gmat(coords_ref1,cutoff).reshape(-1)
    coords_ref2 = ref.xyz[0,nn_ref.indeces_lcs[:,l1:]]
    ref_mat2 = tools.calc_gmat(coords_ref2,cutoff).reshape(-1)
   
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


    lcs_idx = nn_traj.indeces_lcs
    idxs_combo = []
    results = []
    count = 1
    for i in xrange(traj.n_frames):

        # calculate eRMSD for strand1 
        gmats1 = [tools.calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1) for j in all_idx1]
        dd1 = cdist([ref_mat1],gmats1)
        low1 = np.where(dd1[0]<treshold*np.sqrt(l1))
        
        # calculate eRMSD for strand2 
        gmats2 = [tools.calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1) for j in all_idx2]
        dd2 = cdist([ref_mat2],gmats2)
        low2 = np.where(dd2[0]<treshold*np.sqrt(l2))

        # do combination
        combo_list = list(itertools.product(low1[0],low2[0]))
        
        gmats_combo = []
        
        for cc in combo_list:
            llc = len(set(all_idx1[cc[0]] + all_idx2[cc[1]]))
            # skip overlapping
            if(llc != l1 + l2): continue
            
            # skip distant
            com1 = np.average(np.average(traj.xyz[i,lcs_idx[:,all_idx1[cc[0]]]],axis=0),axis=0)
            com2 = np.average(np.average(traj.xyz[i,lcs_idx[:,all_idx2[cc[1]]]],axis=0),axis=0)
            dcoms = np.sqrt(np.sum((com1-com2)**2))
            if(dcoms > 2.5*dcom): continue

            idx_combo = all_idx1[cc[0]] + all_idx2[cc[1]]
            idxs_combo.append(idx_combo)
            gmats_combo.append(tools.calc_gmat(traj.xyz[i,lcs_idx[:,idx_combo]],cutoff).reshape(-1))

        # calculate distances
        dd_combo = cdist([ref_mat],gmats_combo)
        low_combo = np.where(dd_combo[0]<treshold*np.sqrt(l1 + l2))

        for k in low_combo[0]:
            
            #print idxs_combo[k]
            resname_idxs = [nn_traj.rna_seq[l] for l  in idxs_combo[k]]
            
            results.append([i,dd_combo[0,k]/np.sqrt(l1 + l2),resname_idxs])

            #print results[-1]
            # Write aligned PDB 
            if(write != None):
                pdb_out = "%s_%05d_%s_%d.pdb" % (write,count,resname_idxs[k][0],i)
                # slice trajectory
                tmp_atoms = []
                tmp_res =[]
                for r1 in idxs_combo[k]:
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

