#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2017 Sandro Bottaro (sandro.bottaro@bio.ku.dk)
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import, division, print_function

# Make sure that range returns an iterator also in python2 (using future module)
from builtins import range

import sys
import mdtraj as md
from scipy.spatial import distance
import itertools
import numpy as np
import os
from . import definitions
from . import nucleic
from . import calc_mats as ff
    
def ermsd(reference,target,cutoff=2.4,topology=None):
    
    """
    Calculate ermsd between reference and target structures  

    Parameters
    ----------
    reference : string 
         Filename of reference structure, any format accepted by MDtraj can be used.
    target : string 
         Filename of target structure. If a trajectory is provided, a topology file must be specified.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    cutoff :  float, optional
         Cutoff for eRMSD calculation. 
         The default value of 2.4 should work in most cases. This cutoff value roughly correspond to considering pair of bases whose distance is within an ellipsoidal cutoff with axis x=y=2.4*5 = 12 Angstrom and z=2.4*3=7.2 Angstrom. Larger values of cutoff can be useful when analyzing unstructured/flexible molecules.
    Returns
    -------
        array :
            eRMSD distance numpy array with dimension *m*,  the number of structures in target.
    
    """

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
        
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return ermsd_traj(ref,traj,cutoff=cutoff)

def ermsd_traj(reference,traj,cutoff=2.4):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = reference.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))
    
    coords_ref = reference.xyz[0,nn_ref.indeces_lcs]
    ref_mat = ff.calc_gmat(coords_ref,cutoff).reshape(-1)
    #rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    gmats = []
    for i in range(traj.n_frames):
        coords_lcs = traj.xyz[i,nn_traj.indeces_lcs]
        gmats.append(ff.calc_gmat(coords_lcs,cutoff).reshape(-1))
    dd = distance.cdist([ref_mat],gmats)/np.sqrt(len(nn_traj.ok_residues))
    return dd[0]

############## ERMSD ###############


def dump_rvec(filename,topology=None,cutoff=2.4):
    """
    Calculate relative position of pair of nucleobases within ellipsoidal cutoff

    Parameters
    ----------
    reference : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    cutoff :  float, optional
         Cutoff for eRMSD calculation. 
         This cutoff value roughly correspond to considering pair of bases whose distance is within an ellipsoidal cutoff with axis x=y=2.4*5 = 12 Angstrom and z=2.4*3=7.2 Angstrom. Larger values of cutoff can be useful when analyzing unstructured/flexible molecules.
    Returns
    -------
    rmat :
        Numpy array with dimension (m,n,n,3). *m* is the number of structures in target, *n* is the number of nucleotides. As an example, the position of base 10 in the reference system of base 9 in the fourth frame is given by v = rmat[3,8,9], where v is a 3-dimensional vector
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX

    """

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return  dump_rvec_traj(traj,cutoff=cutoff)

def dump_rvec_traj(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    rvecs = []
    for i in range(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        rvecs.append(ff.calc_rmat(coords_lcs,cutoff))
        
    return np.asarray(rvecs), nn.rna_seq

###############################################

def dump_gvec(filename,topology=None,cutoff=2.4):
    
    """
    Calculate relative position of pair of nucleobases within ellipsoidal cutoff

    Parameters
    ----------
    reference : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    cutoff :  float, optional
         Cutoff for eRMSD calculation. 
         This cutoff value roughly correspond to considering pair of bases whose distance is within an ellipsoidal cutoff with axis x=y=2.4*5 = 12 Angstrom and z=2.4*3=7.2 Angstrom. Larger values of cutoff can be useful when analyzing unstructured/flexible molecules.
    Returns
    -------
    gmat :
        Numpy array with dimension (m,n,n,4). *m* is the number of structures in target, *n* is the number of nucleotides. As an example, the position of base 10 in the reference system of base 9 in the fourth frame is given by v = rmat[3,8,9], where v is a 4-dimensional vector
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX
       

    """
    
    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return dump_gvec_traj(traj,cutoff=cutoff)

def dump_gvec_traj(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    gvecs = []
    for i in range(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        gvecs.append(ff.calc_gmat(coords_lcs,cutoff))
    return np.asarray(gvecs), nn.rna_seq

#################################################


def rmsd(reference,target,topology=None,out=None):
    
    """
    Calculate rmsd after optimal alignment between reference and target structures. Superposition and RMSD calculations are performed using all heavy atoms. 
    If the sequence of reference and target is different, only backbone/sugar heavy atoms are used.

    Parameters
    ----------
    reference : string 
         Filename of reference structure, any format accepted by MDtraj can be used.
    target : string 
         Filename of target structure. If a trajectory is provided, a topology file must be specified.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    out :  string, optional
         If a string is specified, superimposed PDB structures are written to disk with the specified prefix.
    Returns
    -------
    array :
        RMSD distance (in nm) numpy array with dimension *m*,  the number of structures in target.
    
    """

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target

    return rmsd_traj(ref,traj,out=out)


def rmsd_traj(reference,traj,out=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = reference.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)
    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))

    # loop over residues and find common heavy atoms
    idx_ref = []
    idx_target = []
    
    for ii in range(len(nn_ref.ok_residues)):
        
        res1 = nn_ref.ok_residues[ii]
        res2 = nn_traj.ok_residues[ii]

        resname1 = nn_ref.rna_seq_id[ii]
        resname2 = nn_traj.rna_seq_id[ii]
        
        # if the nucleotide is the same, use all atoms
        if(res1.name == res2.name):
            
            name2 = [at.name for at in res2.atoms if  at.name in definitions.nt_atoms[resname2]]
            #name1 = [at.name for at in res1.atoms if  at.name in definitions.nb_atoms[res1.name]]
            #print name2,name1
            for at in res1.atoms:
                if at.name in definitions.nt_atoms[resname1]:
                    if(at.name in name2):
                        #print at.name
                        idx_ref.append(at.index)
                        idx_target.append(((res2.atom(at.name)).index))
        # else, use bb only
        else:
            name2 = [at.name for at in res2.atoms if  at.name in definitions.bb_atoms]
            for at in res1.atoms:
                if at.name in definitions.bb_atoms:
                    if(at.name in name2):
                        idx_ref.append(at.index)
                        idx_target.append(((res2.atom(at.name)).index))
            
    print("# found ",len(idx_ref), "atoms in common")
    
    if(len(idx_ref)<3):
        warn =  "# Only  %d atoms in common. abort.\n" % len(idx_ref)
        sys.stderr.write(warn)
        sys.exit(1)
        
    traj.superpose(reference,atom_indices=idx_target, ref_atom_indices=idx_ref)
    if(out!=None):
        traj.save(out)
    rmsd = np.sqrt(3*np.mean((traj.xyz[:, idx_target, :] - reference.xyz[0,idx_ref, :])**2, axis=(1,2)))
    return rmsd

########################################################

def backbone_angles(filename,topology=None,residues=None,angles=None):

    """
    Calculate backbone ([alpha,beta,gamma,delta,espilon,zeta]) and glycosydic (chi) torsion angles.

    Parameters
    ----------
    reference : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    residues :  list, optional
         If a list of residues is specified, only the selected residues will be calculated. Otherwise, the calculation is performed for all residues.
         The residue naming convention is RESNAME_RESNUMBER_CHAININDEX
    angles : list, optional
         If a list of angles is specified, only the selected angles will be calculated. 
         Otherwise, the calculation is performed for all torsion angles. 
    Returns
    -------
    array :
        Torsion angles in radians. A Numpy array with dimension (m,n,q) is returned. *m* is the number of structures in target, *n* is the number of residues, and q is the number of requested angles (7 by default).
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX

    """
    
    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return backbone_angles_traj(traj,residues=residues,angles=angles)

def backbone_angles_traj(traj,residues=None,angles=None):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =  nn.get_bb_torsion_idx(residues)
               
    if(angles==None):
        idx_angles = np.arange(all_idx.shape[1])
    else:
  
        idx_angles = []
        for i in range(len(angles)):
            if(angles[i] in definitions.bb_angles):
                idx_angles.append(definitions.bb_angles.index(angles[i]))
            else:
                msg = "# Fatal error. requested angle \"%s\" not available.\n" % angles[i]
                msg += "# Choose from: %s \n" % definitions.bb_angles
                sys.stderr.write(msg)
                sys.exit(1)
   

    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    
    torsions = md.compute_dihedrals(traj,idxs,opt=True)
    
    # set to NaN where atoms are missing
    torsions[:,np.where(np.sum(idxs,axis=1)==0)[0]] = np.nan
    
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))
    
    return torsions, rr
########################################################
    
def sugar_angles(filename,topology=None,residues=None,angles=None):
    
    """
    Calculate sugar [nu1,nu2,nu3,nu4,nu5] torsion angles.

    Parameters
    ----------
    reference : string  
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional 
         Topology filename. Must be specified if target is a trajectory.
    residues :  list, optional
         If a list of residues is specified, only the selected residues will be calculated. Otherwise, the calculation is performed for all residues.
         The residue naming convention is RESNAME_RESNUMBER_CHAININDEX
    angles : list, optional
         If a list of angles is specified, only the selected angles will be calculated. 
         Otherwise, the calculation is performed for all torsion angles. 

    Returns
    -------
    array :
        Torsion angles in radians. A Numpy array with dimension (m,n,q) is returned. *m* is the number of structures in target, *n* is the number of residues, and q is the number of requested angles (5 by default).
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX

    """

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return sugar_angles_traj(traj,residues=residues,angles=angles)

def sugar_angles_traj(traj,residues=None,angles=None):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =  nn.get_sugar_torsion_idx(residues)
    if(angles==None):
        idx_angles = np.arange(all_idx.shape[1])
    else:
        # find indeces corresponding to angles
        idx_angles = []
        for i in range(len(angles)):
            if(angles[i] in definitions.sugar_angles):
                idx_angles.append(definitions.sugar_angles.index(angles[i]))
            else:
                msg = "# Fatal error. requested angle \"%s\" not available.\n" % angles[i]
                msg += "# Choose from: %s \n" % (definitions.sugar_angles)
                sys.stderr.write(msg)
                sys.exit(1)

    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)

    torsions = md.compute_dihedrals(traj,idxs,opt=True)
    # set to NaN where atoms are missing
    torsions[:,missing[0]] = np.nan
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))

    return torsions, rr

#############################################################

def pucker_angles(filename,topology=None,residues=None):
    
    """
    Calculate sugar pucker pseudorotation  torsion angles: phase and amplitude

    Parameters
    ----------
    reference : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    residues :  list, optional
         If a list of residues is specified, only the selected residues will be calculated. Otherwise, the calculation is performed for all residues.
         The residue naming convention is RESNAME_RESNUMBER_CHAININDEX
    Returns
    -------
    array :
        Phase and amplitude. A Numpy array with dimension (m,n,2) is returned. *m* is the number of structures in target, *n* is the number of residues.
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX

    """

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return pucker_angles_traj(traj,residues=residues)

def pucker_angles_traj(traj,residues=None):

    torsions,rr = sugar_angles_traj(traj,residues=residues)
    x1 = torsions[:,:,4] +  torsions[:,:,1] -  torsions[:,:,3] -   torsions[:,:,0]
    x2 = 3.0776835*torsions[:,:,2]
    phase = np.arctan2(x1,x2)
    phase[np.where(phase<0.0)] += 2.0*np.pi
    tm = torsions[:,:,2]/np.cos(phase)
    angles = np.dstack((phase,tm))
    return angles, rr

################################################################

################################################################

def jcouplings(filename,topology=None,residues=None,couplings=None,raw=False):
    
    """
    Calculate 3J scalar couplings from structure using the Karplus equations.

    Parameters
    ----------
    reference : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    residues :  list, optional
         If a list of residues is specified, only the selected residues will be calculated. Otherwise, the calculation is performed for all residues.
         The residue naming convention is RESNAME_RESNUMBER_CHAININDEX
    couplings:
         If a list of couplings to be chosen from [H1H2, H2H3,H3H4,1H5P,2H5P,C4Pb,1H5H4,2H5H4,H3P,C4Pe,H1C2/4,H1C6/8] 
         is specified, only the selected couplings will be calculated. 
         Otherwise, the calculation is performed for all of them. 
    raw: bool, optional
         raw values of the angles are returned. 

    Returns
    -------
    array :
         Scalar couplings in Hz. A Numpy array with dimension (m,n,q) is returned. *m* is the number of structures in target, *n* is the number of residues and *q* is the number of couplings. If raw = True, a (m,n,5) array is returned with the values of the torsion angles used for the calculation.
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX

    """

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return jcouplings_traj(traj,residues=residues,couplings=couplings,raw=raw)

def jcouplings_traj(traj,residues=None,couplings=None,raw=False):
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    all_idx,rr =nn.get_coupling_idx(residues)

    # if not provided, calculate all couplings
    if(couplings==None):
        idx_angles = np.arange(all_idx.shape[1])
        couplings = [str(el) for el in definitions.couplings_idx.keys()]
        idx_angles1 = [el for el in definitions.couplings_idx.values()]
    else:
        idx_angles = []
        idx_angles1 = []
        for i in range(len(couplings)):
            if(couplings[i] in definitions.couplings_idx.keys()):
                idx_angles.append(definitions.couplings_idx[couplings[i]])
                idx_angles1.append(i)
            else:
                msg = "# Fatal error. requested coupling \"%s\" not available.\n" % couplings[i]
                msg += "# Choose from: %s \n" % couplings_idx
                sys.stderr.write(msg)
                sys.exit(1)
   
    idxs = (all_idx[:,idx_angles,:]).reshape(-1,4)
    missing = np.where(np.sum(idxs,axis=1)==0)
    
    torsions = md.compute_dihedrals(traj,idxs,opt=True)
    
    # set to NaN where atoms are missing
    torsions[:,np.where(np.sum(idxs,axis=1)==0)[0]] = np.nan
    
    torsions = torsions.reshape((traj.n_frames,all_idx.shape[0],len(idx_angles)))

    #jcouplings = np.zeros((torsions.shape[0],torsions.shape[1],len(couplings)))
    # now calculate couplings
    if(raw):
        return torsions,rr

    jcouplings = np.empty((traj.n_frames,all_idx.shape[0],len(couplings)))*np.nan

    for i in range(len(couplings)):
        # get karplus coefficients
        
        coef = definitions.couplings_karplus[couplings[i]]
        #ii = definitions.couplings_idx[couplings[i]]
        ii =  idx_angles1[i]

        #print ii, couplings[i],
        angles = np.copy(torsions[:,:,ii])
        
        # add phase
        if(coef[4] != 0.0):
            angles += coef[4]
        cos = np.cos(angles)
        val = coef[0]*cos*cos + coef[1]*cos
        # add shift
        if(coef[2] != 0.0):
            val += coef[2]
        # add generalized karplus term
        if(coef[3] != 0.0):
            sin = np.sin(angles)
            val += coef[3]*cos*sin
        jcouplings[:,:,i] = val
    
    return jcouplings,rr


##############################################################

def ss_motif(query,target,topology=None,threshold=0.8,cutoff=2.4,sequence=None,out=None,bulges=0):
    
    """
    Find single stranded motif similar to *query* in *target*

    Parameters
    ----------
    query : string 
         Filename of query structure, any format accepted by MDtraj can be used.
    target : string 
         Filename of target structure. If a trajectory is provided, a topology file must be specified.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    threshold : float, optional
         all substructures in target with eRMSD < threshold will be returned. 
    sequence: string, optional
         By default, the search is performed in a sequence-independent manner, unless a specific sequence is specified. Abbreviations (R/Y/N) are accepted.
    out: string, optional
         Hits are written to PDB files and aligned to query with the specified prefix. If *out* is not specified, PDB are not written.
    bulges: int, optional
         Maximum number of allowed bulges, i.e. maximium number of inserted nucleotides. Default value is 0.
    cutoff :  float, optional
         Cutoff for eRMSD calculation. 
         The default value of 2.4 should work in most cases. This cutoff value roughly correspond to considering pair of bases whose distance is within an ellipsoidal cutoff with axis x=y=2.4*5 = 12 Angstrom and z=2.4*3=7.2 Angstrom. Larger values of cutoff can be useful when analyzing unstructured/flexible molecules.
    
    Returns
    -------
        list :
            list of results. Each element in the list contain: the frame number, the eRMSD from query and the residues
    
    """

    ref = md.load(query)
    warn =  "# Loaded query %s \n" % query
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)
    
    return ss_motif_traj(ref,traj,threshold=threshold,cutoff=cutoff,sequence=sequence,out=out,bulges=bulges)

def ss_motif_traj(ref,traj,threshold=0.8,cutoff=2.4,sequence=None,bulges=0,out=None):
    
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
    ref_mat = ff.calc_gmat(coords_ref,cutoff).reshape(-1)
    
    rna_seq = nn_traj.rna_seq_id
    res_idxs = definitions.get_idx(rna_seq,sequence,bulges)
    resname_idxs = [[nn_traj.rna_seq[l] for l in rr]  for rr  in res_idxs]
    if(len(res_idxs)==0):
        return []

    lcs_idx = nn_traj.indeces_lcs
    results = []
    count = 1
    for i in range(traj.n_frames):
        
        gmats = [ ff.calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1)  for j in res_idxs]
            
        dd = distance.cdist([ref_mat],gmats)/np.sqrt(ll)
        low = np.where(dd[0]<threshold)
        for k in low[0]:
            results.append([i,dd[0,k],resname_idxs[k]])

            # Write aligned PDB 
            if(out != None):
                pdb_out = "%s_%05d_%s_%d.pdb" % (out,count,resname_idxs[k][0],i)
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

##########################################################################################

def ds_motif(query,target,l1,l2,threshold=0.9,cutoff=2.4,topology=None,sequence=None,bulges=0,out=None):
    
    """
    Find single stranded motif similar to *query* in *target*

    Parameters
    ----------
    query : string 
         Filename of query structure, any format accepted by MDtraj can be used.
    target : string 
         Filename of target structure. If a trajectory is provided, a topology file must be specified.
    l1 : int 
         Number of nucleotides in the first strand
    l2 : int 
         Number of nucleotides in the second strand
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    threshold : float, optional
         all substructures in target with eRMSD < threshold will be returned. 
    sequence: string, optional
         By default, the search is performed in a sequence-independent manner, unless a specific sequence is specified. Abbreviations (R/Y/N) are accepted.
    out: string, optional
         Hits are written to PDB files and aligned to query with the specified prefix. If *out* is not specified, PDB are not written.
    bulges: int, optional
         Maximum number of allowed bulges, i.e. maximium number of inserted nucleotides. Default value is 0.
    cutoff :  float, optional
         Cutoff for eRMSD calculation. 
         The default value of 2.4 should work in most cases. This cutoff value roughly correspond to considering pair of bases whose distance is within an ellipsoidal cutoff with axis x=y=2.4*5 = 12 Angstrom and z=2.4*3=7.2 Angstrom. Larger values of cutoff can be useful when analyzing unstructured/flexible molecules.
    
    Returns
    -------
        list :
            list of results. Each element in the list contain: the frame number, the eRMSD from query and the residues
    
    """

    ref = md.load(query)
    warn =  "# Loaded query %s \n" % query        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return ds_motif_traj(ref,traj,l1,l2,threshold=threshold,sequence=sequence,bulges=bulges,out=out)

def ds_motif_traj(ref,traj,l1,l2,threshold=0.9,cutoff=2.4,sequence=None,bulges=0,out=None):
    
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
    ref_mat = ff.calc_gmat(coords_ref,cutoff).reshape(-1)

    coords_ref1 = ref.xyz[0,nn_ref.indeces_lcs[:,0:l1]]

    ref_mat1 = ff.calc_gmat(coords_ref1,cutoff).reshape(-1)
    coords_ref2 = ref.xyz[0,nn_ref.indeces_lcs[:,l1:]]
    ref_mat2 = ff.calc_gmat(coords_ref2,cutoff).reshape(-1)
   
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
    for i in range(traj.n_frames):

        # calculate eRMSD for strand1 
        gmats1 = [ff.calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1) for j in all_idx1]
        dd1 = distance.cdist([ref_mat1],gmats1)
        low1 = np.where(dd1[0]<threshold*np.sqrt(l1))
        
        # calculate eRMSD for strand2 
        gmats2 = [ff.calc_gmat(traj.xyz[i,lcs_idx[:,j]],cutoff).reshape(-1) for j in all_idx2]
        dd2 = distance.cdist([ref_mat2],gmats2)
        low2 = np.where(dd2[0]<threshold*np.sqrt(l2))

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
            gmats_combo.append(ff.calc_gmat(traj.xyz[i,lcs_idx[:,idx_combo]],cutoff).reshape(-1))

        # calculate distances
        dd_combo = distance.cdist([ref_mat],gmats_combo)
        low_combo = np.where(dd_combo[0]<threshold*np.sqrt(l1 + l2))

        for k in low_combo[0]:
            
            #print idxs_combo[k]
            resname_idxs = [nn_traj.rna_seq[l] for l  in idxs_combo[k]]
            
            results.append([i,dd_combo[0,k]/np.sqrt(l1 + l2),resname_idxs])

            #print results[-1]
            # Write aligned PDB 
            if(out != None):
                pdb_out = "%s_%05d_%s_%d.pdb" % (out,count,resname_idxs[k][0],i)
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


######################################################################################



def annotate(filename,topology=None):
    
    """
    Find base-pair and base-stacking 

    Parameters
    ----------
    filename : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    Returns
    -------
    stackings : list
    pairings : list
        stackings and pairings contains the list of interactions for the N frames in the PDB/trajectory file and it is organized in the following way: for a given frame i=1..n there are k=1..q interactions between residues with index pairings[i][0][k][0] and pairings[i][0][k][1]. The type of interaction is specified at the element pairings[i][1][k].
    seq : 
        List of residue names. Each residue is identified with the string RESNAME_RESNUMBER_CHAININDEX       

    """
    
    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    
    return annotate_traj(traj)


def annotate_traj(traj):

    # this is the binning for annotation
    bins = [0,1.84,3.84,2.*np.pi]
    bins_label = ["W","H","S"]
    
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    
    #max_r  = np.max(definitions.f_factors)*1.58
    #condensed_idx =  np.triu_indices(len(nn.ok_residues), 1)

    stackings = []
    pairings = []
    
    for i in range(traj.n_frames):

        # calculate LCS
        coords = traj.xyz[i,nn.indeces_lcs]

        # find bases in close contact (within ellipsoid w radius 1.7)
        pairs,vectors,angles = ff.calc_mat_annotation(coords)
        if(len(pairs)==0):
            stackings.append([[],[]])
            pairings.append([[],[]])
            continue
        
        # calculate rho
        rho_12 = vectors[:,0,0]**2 + vectors[:,0,1]**2
        rho_21 = vectors[:,1,0]**2 + vectors[:,1,1]**2

        # calculate z squared 
        z_12 = vectors[:,0,2]**2
        z_21 = vectors[:,1,2]**2
        # find stacked bases 

        # z_ij  AND z_ji > 2 AA
        stackz_12 = np.where(z_12>0.04)
        stackz_21 = np.where(z_21>0.04)
        
        # rho_ij OR rho_ji < 2.5 AA
        rhoz_12 =  np.where(rho_12<0.0625)
        rhoz_21 =  np.where(rho_21<0.0625)
        union_rho = np.union1d(rhoz_12[0],rhoz_21[0])

        # angle between normal planes < 40 deg
        anglez = np.where(np.abs(angles) >0.766)

        # intersect all criteria
        inter_zeta = np.intersect1d(stackz_12,stackz_21)
        inter_stack = np.intersect1d(union_rho,inter_zeta)
        stacked = np.intersect1d(inter_stack,anglez)

        # now I have found all stackings.
        #stacked_pairs = [[nn.rna_seq[pairs[k,0]],nn.rna_seq[pairs[k,1]]] for k in stacked]
        stacked_pairs = [[pairs[k,0],pairs[k,1]] for k in stacked]
        stacked_annotation = np.chararray((len(stacked_pairs),2), unicode='True')
        stacked_annotation[:] = ">"
        
        # revert where z_ij is negative
        rev1 = np.where(vectors[stacked,0,2]<0)
        stacked_annotation[rev1[0],0] = "<"

        rev2 = np.where(vectors[stacked,1,2]>0)
        stacked_annotation[rev2[0],1] = "<"

        stacked_annotation = ["".join(el) for el in stacked_annotation]
        ######################################################

        # find paired bases  (z_ij < 2 AA AND z_j1 < 2 AA)
        #paired = [j for j in range(len(pairs)) if((j not in stackz_12[0]) and (j not in stackz_21[0]))]
        paired = [j for j in range(len(pairs)) if((j not in stackz_12[0]) or (j not in stackz_21[0]))]
        #paired_pairs = [[nn.rna_seq[pairs[k,0]],nn.rna_seq[pairs[k,1]]] for k in paired]
        paired_pairs = [[pairs[k,0],pairs[k,1]] for k in paired]
        
        # calculate edge angle. subtract 0.16 as Watson edge is not zero
        edge_angles_1 = np.arctan2(vectors[paired,0,1],vectors[paired,0,0]) - definitions.theta1
        edge_angles_2 = np.arctan2(vectors[paired,1,1],vectors[paired,1,0]) - definitions.theta1

        # shift to 0-2pi range
        edge_angles_1[np.where(edge_angles_1<0.0)] += 2.*np.pi
        edge_angles_2[np.where(edge_angles_2<0.0)] += 2.*np.pi

        # find edge: 0 Watson, 1:Hoogsteen, 2sugar
        edge_1 = np.digitize(edge_angles_1,bins)-1
        edge_2 = np.digitize(edge_angles_2,bins)-1

        
        paired_annotation = np.chararray((len(paired_pairs),3), unicode=True)
        paired_annotation[:] = "X"
        
        # now explicit loop, a bit messy. sorry, Guido.
        for j in range(len(paired_pairs)):

            # if angle is larger than 60 deg, skip
            if(np.abs(angles[paired[j]]) < 0.5 ): continue

            index1 = paired_pairs[j][0]
            index2 = paired_pairs[j][1]
            
            # find donor and acceptor in first base
            r1_donor = nn.donors[index1]
            r1_acceptor = nn.acceptors[index1]
            # find donor and acceptor in second base
            r2_donor = nn.donors[index2]
            r2_acceptor = nn.acceptors[index2]
            combo_list = list(itertools.product(r1_donor,r2_acceptor)) + list(itertools.product(r1_acceptor,r2_donor))
            
            # distances between donor and acceptors
            delta = np.diff(traj.xyz[i,combo_list],axis=1)
            dist_sq = np.sum(delta**2,axis=2)
            # number or distances less than 3.3 AA is n_hbonds
            n_hbonds = (dist_sq<0.1089).sum()
            # if no hydrogen bonds, skip 
            if(n_hbonds==0): continue

            # find edge
            paired_annotation[j][0] = bins_label[edge_1[j]]
            paired_annotation[j][1] = bins_label[edge_2[j]]

            gidxs = [nn.indeces_glyco[index1][1],nn.indeces_glyco[index1][0],nn.indeces_glyco[index2][0],nn.indeces_glyco[index2][1]]

            # if atoms in glyco are missing, do not calculate cis/trans
            if(None in gidxs):
                paired_annotation[j][2] = "x"
            else:
                angle_glyco = ff.dihedral(traj.xyz[i,gidxs[0]],traj.xyz[i,gidxs[1]],traj.xyz[i,gidxs[2]],traj.xyz[i,gidxs[3]])
                if(np.abs(angle_glyco) > 0.5*np.pi):
                    paired_annotation[j][2] = "t"
                else:
                    paired_annotation[j][2] = "c"
                
            # if is WWc, check for Watson-crick and GU
            if("".join(paired_annotation[j]) == "WWc"):
                r1 =  nn.rna_seq_id[index1]
                r2 =  nn.rna_seq_id[index2]
                ll = "".join(sorted([r1,r2]))
                if((z_12[paired[j]] < 0.04) and (z_21[paired[j]]<0.04)):
                    if(((ll=="AU" and n_hbonds > 1) or ( ll == "CG" and n_hbonds > 2))):
                        paired_annotation[j] = ["W","C","c"]
                    if(ll=="GU" and n_hbonds > 1):
                        paired_annotation[j] = ["G","U","c"]
                    
        paired_annotation = ["".join(el) for el in paired_annotation]

        
        stackings.append([stacked_pairs,stacked_annotation])
        pairings.append([paired_pairs,paired_annotation])

    return stackings, pairings, nn.rna_seq

## DOT-BRACKET ##

def dot_bracket(pairings,sequence):

    """
    calculate dot-bracket annotation 

    Parameters
    ----------
    pairings : list 
         List of pairings as returned by annotate function
    sequence : list
         Sequence as returned by annotate function
    Returns
    -------
    dot_bracket: list
        strings with the dot-bracket annotation, one per frame
    ss: string
        string with sequence in dbn format

    """
    def get_seq(res):
        ss = ""
        chain_break = []
        for ii,r in enumerate(res):
            el = r.split("_")[0]
            if(el in definitions.residue_dict):
                name = definitions.residue_dict[el][-1]
            else:
                name = "X"
            ss += name
            if(ii!=len(res)-1):
                if(res[ii+1].split("_")[2] != res[ii].split("_")[2]):
                    ss+= "&"
                    chain_break.append(ii)
        ss += "\n"
        return ss,chain_break

    ss,chain_break = get_seq(sequence)
    # loop over frames
    dot_bracket = []
    for k,pp in enumerate(pairings):
        
        # loop over annotation
        string = ["."]*len(sequence)
        for e in range(len(pp[0])):
            if(pp[1][e]=="WCc"):
                idx1 = pp[0][e][0]
                idx2 = pp[0][e][1]

                # check that bracket is not there already. This might happen in simulations
                if(string[idx1]!="." or string[idx2]!="." ):
                    warn = "# Frame %d, residue %s has a double WC base-pairs. " % (k,sequence[idx1],sequence[idx2])
                    warn += " Dot-bracket annotation set to XXX \n"
                    sys.stderr.write(warn)
                    dot_bracket.append("XXX")                            
                    continue
                # decide brackets
                level = 0                    
                for br in definitions.cl:
                    if(string[idx1:idx2].count(br) != 0):
                        level += 1
                    else:
                        break
                string[idx1] = definitions.op[level]
                string[idx2] = definitions.cl[level]
                
        # add chain breaks
        new_string = ""
        for ii,s in enumerate(string):
            new_string += s
            if(ii in chain_break):  new_string += "&"
        dot_bracket.append(new_string)
    #print(dot_bracket)
    #exit()
    return dot_bracket, ss


#############################################################

def snippet(pdb,sequence,outdir=None):

    """
    Extract fragments fro pdb with a givn sequence

    Parameters
    ----------
    pdb : string 
         name of PDB file. Only PDB are accepted. 
    sequence: 
         sequence to extract. N/R/Y abbreviations are accepted. 
    outdir : string
        specify output directory

    """

    from . import reader

    if(outdir==None):
        outdir = os.getcwd()

    atoms = ["C2","C4","C6"]
    
    # check query sequence
    for item in sequence:
        if(item not in definitions.known_abbrev):
            print("# FATAL Error. Symbol ", item, " not known. Use AUCG/NYRSWKMBDHV")
            return 1
        if(item == "%"):
            print("# Fatal error. Single strand only")
            return 1

    ll = [len(el) for el in sequence]
    cur_pdb = reader.Pdb(pdb,res_mode="R",permissive=True)
    cur_len = len(cur_pdb.model.sequence)
    indeces = definitions.get_idx(cur_pdb.model.sequence,sequence,bulges=0)
   
         
    idx = 0
    ii = 0
    while(idx>=0):
      
        for index in indeces:

            # do checks
            skip = False
            for k,res in enumerate(index):
                rr =  cur_pdb.model[res]
                # check that atoms in the base are in place
                out = [rr.get_idx(aa) for aa in atoms]
                if(np.isnan(np.sum(out))):
                    skip = True

                # check connectivity
                if(k<len(index)-1):
                    # check that O3' and P are connected
                    rrp = cur_pdb.model[index[k+1]]
                    dd = np.sqrt(np.sum(( np.array(rrp["P"]) -np.array(rr["O3'"]))**2))
                    if(dd>1.7):
                        skip = True
                if(skip):
                    continue
            
            name_pref = pdb[0:-4].split("/")[-1]
            new_pdb = "%s/%s_%s_%05d.pdb" % (outdir,name_pref,cur_pdb.model.sequence_id[index[0]],ii)
            sys.stderr.write("# Writing PDB %s\n" % new_pdb.split("/")[-1])

            fh_pdb = open(new_pdb,'w')
            fh_pdb.write(cur_pdb.model.string_pdb(index,noP=True,center=True))
            fh_pdb.close()
         
            ii += 1
        idx = cur_pdb.read()
      
    
