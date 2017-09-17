import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist
import tools
import sys
import nucleic
import definitions


def ssmotif_traj(ref,traj,treshold=0.9,cutoff=2.4,sequence=None,bulges=0,write=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = ref.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    ll = len(nn_ref.ok_residues)
    if(sequence==None):
        sequence = "N"*ll
    else:
        assert(len(sequence)==ll)
        
    coords_ref = ref.xyz[0,nn_ref.indeces_lcs]
    ref_mat = tools.calc_gmat(coords_ref,cutoff).reshape(-1)
    
    rna_seq = nn_traj.rna_seq_id
    res_idxs = definitions.get_idx(rna_seq,sequence,bulges)
    resname_idxs = [[nn_traj.rna_seq[l] for l in rr]  for rr  in res_idxs]
    if(len(res_idxs)==0):
        return []

    lcs_idx = nn_traj.indeces_lcs
    results = []
    count = 1
    for i in xrange(traj.n_frames):
        gmats = []
        for j in res_idxs:
            coords_lcs = traj.xyz[i,lcs_idx[:,j]]
            gmats.append(tools.calc_gmat(coords_lcs,cutoff).reshape(-1))
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


def ssmotif(reference,target,treshold=0.9,cutoff=2.4,topology=None,sequence=None,bulges=0,write=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return ssmotif_traj(ref,traj,treshold=treshold,sequence=sequence,bulges=bulges,write=write)

