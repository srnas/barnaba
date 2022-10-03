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
    
def ermsd(reference,target,cutoff=2.4,topology=None,residues_ref=None,residues_target=None):
    
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

    return ermsd_traj(ref,traj,cutoff=cutoff,residues_ref=residues_ref,residues_target=residues_target)

def ermsd_traj(reference,traj,cutoff=2.4,residues_ref=None,residues_target=None):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = reference.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    
    ref_indeces = nn_ref.indeces_lcs
    if(residues_ref!=None): ref_indeces = ref_indeces[:,residues_ref]
        
    coords_ref = reference.xyz[0,ref_indeces]
    ref_mat = ff.calc_gmat(coords_ref,cutoff).reshape(-1)

    #rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    gmats = []
    traj_indeces = nn_traj.indeces_lcs
    if(residues_target!=None): traj_indeces = traj_indeces[:,residues_target]
    #print(residues_target)
    assert ref_indeces.shape[1]==traj_indeces.shape[1], "reference %s, target %s" % (ref_indeces.shape,traj_indeces.shape)
    
    for i in range(traj.n_frames):
        coords_lcs = traj.xyz[i,traj_indeces]
        gmats.append(ff.calc_gmat(coords_lcs,cutoff).reshape(-1))
    dd = distance.cdist([ref_mat],gmats)/np.sqrt(traj_indeces.shape[1])

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


def rmsd(reference,target,topology=None,out=None,heavy_atom=False):
    
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
    heavy_atom :   bool, optional
         If True, all heavy atoms are used for superposition. Default is False

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

    return rmsd_traj(ref,traj,out=out,heavy_atom=heavy_atom)


def rmsd_traj(reference,traj,out=None,heavy_atom=False):
    
    top_traj = traj.topology
    # initialize nucleic class
    nn_traj = nucleic.Nucleic(top_traj)

    top_ref = reference.topology
    # initialize nucleic class
    nn_ref = nucleic.Nucleic(top_ref)

    assert(len(nn_traj.ok_residues)==len(nn_ref.ok_residues))    
    # check that sequence is identical when using heavy atoms
    if(nn_traj.rna_seq_id!=nn_ref.rna_seq_id and heavy_atom==True):
        sys.stderr.write("# Sequences are not identical, cannot superimpose using all heavy atoms. ")
        sys.exit(1)
        
    # loop over residues and find common heavy atoms
    idx_ref = []
    idx_target = []
    
    for ii in range(len(nn_ref.ok_residues)):
        
        res1 = nn_ref.ok_residues[ii]
        res2 = nn_traj.ok_residues[ii]

        resname1 = nn_ref.rna_seq_id[ii]
        resname2 = nn_traj.rna_seq_id[ii]
        
        # if heavy_atom is true, use all atoms
        if(heavy_atom):
            
            name2 = [at.name for at in res2.atoms if  at.name in definitions.nt_atoms[resname2]]
            for at in res1.atoms:
                if at.name in definitions.nt_atoms[resname1]:
                    if(at.name in name2):
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
    if(out!=None): traj.save(out)
    
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

def pucker_angles(filename,topology=None,residues=None,altona=False):
    
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
    if(altona):
        return pucker_altona_traj(traj,residues=residues)
    else:
        return pucker_rao_traj(traj,residues=residues)
        
def pucker_altona_traj(traj,residues=None):

    torsions,rr = sugar_angles_traj(traj,residues=residues)
    x1 = torsions[:,:,4] +  torsions[:,:,1] -  torsions[:,:,3] -   torsions[:,:,0]
    x2 = 3.0776835*torsions[:,:,2]
    phase = np.arctan2(x1,x2)
    phase[np.where(phase<0.0)] += 2.0*np.pi
    tm = torsions[:,:,2]/np.cos(phase)
    angles = np.dstack((phase,tm))
    return angles, rr

def pucker_rao_traj(traj,residues=None):

    torsions,rr = sugar_angles_traj(traj,residues=residues)

    period = 4.*np.pi/5
    A = (2./5.)*(torsions[:,:,0]+np.cos(period)*torsions[:,:,1] + np.cos(2.*period)*torsions[:,:,2] + \
                np.cos(3.*period)*torsions[:,:,3] + np.cos(4.*period)*torsions[:,:,4])
    B = -(2./5.)*(np.sin(period)*torsions[:,:,1] + np.sin(2.*period)*torsions[:,:,2] + \
                np.sin(3.*period)*torsions[:,:,3] + np.sin(4.*period)*torsions[:,:,4])
    
    phase = np.arctan2(B,A) - 0.5*period
    phase[np.where(phase<0.0)] += 2.0*np.pi
    tm = np.sqrt(A*A + B*B)
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
                msg += "# Choose from: %s \n" % definitions.couplings_idx
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
                #for res1,res2 in zip(nn_ref.ok_residues,traj_slice.topology.residues):
                for ii in range(len(nn_ref.ok_residues)):
                    res1 = nn_ref.ok_residues[ii]
                    res2 = list(traj_slice.topology.residues)[ii]
                    
                    name2 = [at.name for at in res2.atoms if  at.name in definitions.bb_atoms]
                    for at in res1.atoms:
                        if at.name in definitions.bb_atoms:
                            if(at.name in name2):
                                idx_ref.append(at.index)
                                idx_target.append(((res2.atom(at.name)).index))
                                #print(at,res2.atom(at.name))

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

    return ds_motif_traj(ref,traj,l1,l2,cutoff=cutoff,threshold=threshold,sequence=sequence,bulges=bulges,out=out)

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

                pdb_out = "%s_%05d_%s_%d.pdb" % (out,count,resname_idxs[0],i)
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
                for ii in range(len(nn_ref.ok_residues)):
                    res1 = nn_ref.ok_residues[ii]
                    res2 = list(traj_slice.topology.residues)[ii]
                    
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



def annotate(filename,topology=None, stacking_rho_cutoff=2.5, stacking_angle_cutoff=40, pairing_angle_cutoff=60):

    """
    Find base-pair and base-stacking 

    Parameters
    ----------
    filename : string 
         Filename of structure, any format accepted by MDtraj can be used.
    topology : string, optional
         Topology filename. Must be specified if target is a trajectory.
    stacking_rho_cutoff : float
         Cutoff for base stacking. Rho distance. Default is 2.5 (unit in angstroms).
    stacking_angle_cutoff : float
         Cutoff for base stacking. Angle between the vectors normal to the planes of the two bases. Default is 40 (unit in degrees).
    stacking_angle_cutoff : float
         Cutoff for base pairing. Angle between the vectors normal to the planes of the two bases. Default is 60 (unit in degrees).
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
    
    return annotate_traj(traj, stacking_rho_cutoff, stacking_angle_cutoff, pairing_angle_cutoff)


def annotate_traj(traj, stacking_rho_cutoff=2.5, stacking_angle_cutoff=40, pairing_angle_cutoff=60):
 
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
        stackz_12 = np.where(z_12>0.04)  # nm^2
        stackz_21 = np.where(z_21>0.04)  # nm^2
        
        # rho_ij OR rho_ji < 2.5 AA
        #rhoz_12 =  np.where(rho_12<0.0625)
        #rhoz_21 =  np.where(rho_21<0.0625)
        rhoz_12 =  np.where(rho_12<stacking_rho_cutoff * stacking_rho_cutoff * 0.1 * 0.1)  # nm^2
        rhoz_21 =  np.where(rho_21<stacking_rho_cutoff * stacking_rho_cutoff * 0.1 * 0.1)  # nm^2
        union_rho = np.union1d(rhoz_12[0],rhoz_21[0])

        # angle between normal planes < 40 deg
        #anglez = np.where(np.abs(angles) >0.766)
        anglez = np.where(np.abs(angles) > np.cos(np.radians(stacking_angle_cutoff)))   # cosine of angle beween the normal vectors constructed on base i and base j

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
            #if(np.abs(angles[paired[j]]) < 0.5 ): continue
            if(np.abs(angles[paired[j]]) < np.cos(np.radians(pairing_angle_cutoff))): continue

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
                if(string[idx1]!="."):
                    warn = "# Frame %d, residue %s has a double WC base-pairs. " % (k,sequence[idx1])
                    warn += " Dot-bracket annotation set to XXX \n"
                    sys.stderr.write(warn)
                    dot_bracket.append("XXX")                            
                    continue
                if(string[idx2]!="." ):
                    warn = "# Frame %d, residue %s has a double WC base-pairs. " % (k,sequence[idx2])
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


############# SEC_STRUCTURE ########################################

def parse_dotbr(dotbra):
    res_bp = []
    basepairs = []
    for k, i in enumerate(dotbra):
        if i in [")", "]", "}"]:
            if i == ")":
                search = "("
            elif i == "]":
                search = "["
            else:
                search = "{"

            res_bp.append(k)
            for k2 in range(k-1)[::-1]:
                j = dotbra[k2]
                if j == search and k2 not in res_bp:
                    basepairs.append([k2, k])
                    res_bp.append(k2)
                    break
    return basepairs            
    
def parse_dotbracket(threshold, file, weights):
    import barnaba.sec_str_constants as secon
    import re 
    print("Parsing file ", file)    
    with open(file) as f:
        f_dotbracket = f.readlines()
        f.close()
    pdbs = []
    ann_lists = []
    ann_list = {}
    res_bp = []
    n_frames = 0
    n_pdbs = 0
    traj = False
    list_base_pairs = []
    sequence = []
    mult_chain = False
    for line in f_dotbracket:
        if line.startswith("#") and line.split()[1] == "sequence":
            sequence = np.array([r.split("_") for r in line.split()[2].split("-")])
            if len(np.unique(sequence[:,2])) > 1:
                print("Found multiple chains, drawing first chain")
          #      sequence = sequence[np.where(sequence[:,2]==np.unique(sequence[:,2])[0])]
          #      print(sequence)
          #      mult_chain = True
        if not line.startswith("#"):
            l = line.split()
            if "pdb" in l[0].lower():
                dotbr = l[1]
                if len(dotbr) != len(sequence):
                    print(mult_chain)
                    if not mult_chain:
                        sys.exit("Dot-bracket length does not match sequence length")
                    else:
                        dotbr = dotbr[:len(sequence)]
                base_pairs = parse_dotbr(dotbr)
                for bp in base_pairs:
                    ann_list[bp[0], bp[1], "WCc"] = 1.
                list_base_pairs = base_pairs    
                ann_lists.append(dict(ann_list))
                ann_list = {}
            else:
                try: int(l[0])
                except ValueError:
                    continue
                else:
                    traj = True
                    dotbr = l[1]
                    if len(dotbr) != len(sequence):
                        if not mult_chain:
                            sys.exit("Dot-bracket length does not match sequence length")    
                        else:
                            dotbr = dotbr[:len(sequence)]
                    n_frames += 1
                    base_pairs = parse_dotbr(dotbr)
                    for bp in base_pairs:
                        if bp not in list_base_pairs:
                            list_base_pairs.append(bp)
                        try:
                            ann_list[bp[0], bp[1], "WCc"]
                        except:
                            if len(weights) == 0:
                                ann_list[bp[0], bp[1], "WCc"] = 1.
                            else:
                                ann_list[bp[0], bp[1], "WCc"] = weights[n_frames-1]
                        else: 
                            if len(weights) == 0:
                                ann_list[bp[0], bp[1], "WCc"] += 1.
                            else:    
                                ann_list[bp[0], bp[1], "WCc"] += weights[n_frames-1]
    if len(list_base_pairs) == 0:
        sys.exit("No basepairs found.")
    if traj:
        for ann, value in ann_list.items():
            if len(weights) == 0:
                ann_list[ann] /= n_frames
            else:
                ann_list[ann] /= sum(weights)
    chains = [0]
    return sequence, ann_list, list_base_pairs, n_frames    

def parse_annotations(threshold, file, weights):
    import barnaba.sec_str_constants as secon
    import re 
    print("Parsing file ", file)    
    with open(file) as f:
        annotation = f.readlines()
        f.close()
    
    re_pdb = re.compile("^\s*#\s+PDB\s+.*$")
    re_frame = re.compile("^\s*#\s+Frame.*$")
    re_ann = re.compile("^\s*([AGCU])_([0-9]+)_([0-9]+)\s+([AGCU])_([0-9]+)_([0-9]+)\s+.*$")
    chains = []
    old_C = "X"
    ann_lists = []
    ann_list = {}
    n_frames = 0
    for line in annotation:
        if line.startswith("#") and line.split()[1] == "sequence":
            sequence = np.array([r.split("_") for r in line.split()[2].split("-")])
        pdb = re_pdb.match(line)
        if pdb:
            n_frames += 1
            continue
        frame = re_frame.match(line)
        if frame:
            n_frames += 1
            continue
        ann = re_ann.match(line)
        if ann:
            cols = line.split()
            i_N = cols[0].split("_")[0]
            ri = int(cols[0].split("_")[1])
            i_C = cols[0].split("_")[2]
            j_N = cols[1].split("_")[0]
            rj = int(cols[1].split("_")[1])
            j_C = cols[1].split("_")[2]
            ann_ij = cols[2]
            if ann_ij == "XXX":
                continue
            i = np.where(np.array(sequence)[:,1].astype(int) == ri)[0][0]
            j = np.where(np.array(sequence)[:,1].astype(int) == rj)[0][0]
            if old_C != "X" and (i_C != j_C or i_C != old_C):
                ann_lists.append(ann_list)
                ann_list = {}
            if i_C not in chains:
                chains.append(i_C)
            old_C = i_C
            try:
                ann_list[i, j, ann_ij]
            except:
                if len(weights) == 0:
                    ann_list[i, j, ann_ij] = 1.
                else:    
                    ann_list[i, j, ann_ij] = weights[n_frames-1]
            else: 
                if len(weights) == 0:
                    ann_list[i, j, ann_ij] += 1.
                else:    
                    ann_list[i, j, ann_ij] += weights[n_frames-1]
        
    pairs = []
    for ann, value in ann_list.items():
        if len(weights) == 0:
            ann_list[ann] /= n_frames
        else:
            ann_list[ann] /= sum(weights)
        if ann_list[ann] >= threshold and [ann[0], ann[1]] not in pairs:
            pairs.append([ann[0], ann[1]])
    return sequence, ann_list, pairs, n_frames
    
def stems(param_wc, param_bp, param_stack):
    
    #pairs_wc =  np.unique(np.array(param_wc)[:,1:3], axis=0)
    #pairs_all = np.unique(np.array(param_bp + param_stack + param_wc)[:,1:3], axis=0)
    pairs_wc = definitions.unique_rows(np.array(param_wc)[:,1:3])
    pairs_all = definitions.unique_rows(np.array(param_bp + param_stack + param_wc)[:,1:3])
    sums_all = np.unique( np.sum(pairs_all, axis=1))
    ds_all = []
   
    for p in pairs_all:
        added = False
        for di, ds in enumerate(ds_all):
            if abs(ds[-1,0]-p[0])==1. and abs(ds[-1,1]-p[1])==1.:
                if len(ds) > 1:
                    if ds[-1,0]-p[0] == ds[-2,0]-ds[-1,0] and ds[-1,1]-p[1] == ds[-2,1]-ds[-1,1]:
                        added = True
                        ds_all[di] = np.append(ds, [p], axis=0)
                        break
                else:
                    ds_all[di] = np.append(ds, [p], axis=0)
                    added = True
                    break
        if not added:
            ds_all.append(np.array([p]))

    wc_per_ds = [[p for p in s if tuple(p) in set(tuple(x) for x in pairs_wc)] for s in ds_all]
    # sort double strands first by number of wc pairs, then number of all pairs
    sortkeys = np.lexsort((np.array([len(st) for st in ds_all]), np.array([len(wc) for wc in wc_per_ds])))
    ds_all = np.array(ds_all)[sortkeys][::-1]
    stems = []
    diagonal_pairs = []
    lonely_pairs = np.array([])
    for k, d in enumerate(ds_all):
        limits = np.array([[d[:,0].min(), d[:,1].max()], [d[:,0].max(), d[:,1].min()]])
        add = True
        for d2 in ds_all[:k]:
            l2 = np.array([[d2[:,0].min(), d2[:,1].max()], [d2[:,0].max(), d2[:,1].min()]])
            if ((l2[0][0] <= limits[0][0] <= l2[1][0] and l2[1][1] <= limits[0][1] <= l2[0][1]) or
                    (l2[0][0] <= limits[1][0] <= l2[1][0] and l2[1][1] <= limits[1][1] <= l2[0][1])):
                add = False
                diagonal_pairs += list(d)        
        if add:
            if len(d) > 1:
                stems.append(d)
            else:
                if len(lonely_pairs) == 0:
                    lonely_pairs = np.array(d)
                else:    
                    lonely_pairs = np.append(lonely_pairs, np.array(d), axis=0)
    return stems, diagonal_pairs, lonely_pairs

def parameters(pairs, ann_list, n, threshold):
    import barnaba.sec_str_constants as secon

    # parameters:
    # potential type 0: harmonic potential
    # potential type 1: semiharmonic potential to reject bases
    # potential type 2: angular potential to keep ds regions in order
    # potential type 3: semiharmonic potential to reject terminal bases in x direction
  

    param_seq = np.empty((0, 5))
    for i in range(n-1):
        param_seq = np.append(param_seq, [[0, i, i+1, secon.k_seq, secon.d_seq]], axis=0)
    param_bp = []
    param_wc = []
    param_stack = []
    for ann, value in ann_list.items():
        i = ann[0]
        j = ann[1]
        ann_ij = ann[2]
        if value < threshold or abs(i-j) < 2:
            continue
        if ann_ij in secon.list_bp_ct:
            if ann_ij in secon.list_wc_pairs: 
                param_wc.append([0, i, j, secon.k_wc*value, secon.d_short])
            else:
                param_bp.append([0, i, j, secon.k_bp*value, secon.d_short])
        if ann_ij in secon.list_stackings:
            if abs(i-j) == 2:    
                param_stack.append([0, i, j, secon.k_seq*value, secon.d_stack])
            else:
                param_stack.append([0, i, j, secon.k_stack*value, secon.d_short])
    pairs_stems, diagonal_pairs, lonely_pairs = stems(param_wc, param_bp, param_stack)
    param_wc_ds = []
 #   print("pairs_stems", pairs_stems)
 #   print("diagonal_pairs", diagonal_pairs)
 #   print("lonely_pairs", lonely_pairs)
    # Longest stem vertical
    param_stem = np.empty((0, 4))
    if len(pairs_stems) > 0:
        diff = np.absolute(pairs_stems[0][:,1]-pairs_stems[0][:,0])
        keys = np.argsort(diff)
        for k in keys:
            pair = pairs_stems[0][k]
            if abs(pair[0]-pair[1]) > 2:
                param_stem = np.append(param_stem, [[4, int(pair[0:2].min()), int(pair[0:2].max()), 0]], axis=0)
    else:
        if len(lonely_pairs) > 0:
            for pair in lonely_pairs:
                if sum(pair) == sum(lonely_pairs[np.argsort(lonely_pairs[:,0])[0]]):
                    param_stem = np.append(param_stem, [[4, int(pair[0:2].min()), int(pair[0:2].max()), 0]], axis=0)

    for pair in diagonal_pairs:
        for pi in param_wc + param_bp + param_stack:
            if tuple(pi[1:3]) in [tuple(pair), tuple(pair[::-1])]:
                pi[3] = 0
                pi[4] = secon.d_long

    param_ds = []
    param_ang = np.empty((0, 7))
    if len(lonely_pairs) > 0:
        all_stems = pairs_stems
        for p in lonely_pairs:
            all_stems.append(np.array([p]))
    else:
        all_stems = pairs_stems
    for k, stem in enumerate(all_stems):
        stem_limits = np.array([stem[:,0].min(), stem[:,1].max(), stem[:,0].max(), stem[:,1].min()]).astype(int)
        same_stem = stem.copy()
        for k2, s in enumerate(all_stems):
            if k!=k2:
                if sum(s[0]) == sum(stem[0]):
                    same_stem = np.append(same_stem, s, axis=0)
                elif abs(sum(s[0]) - sum(stem[0])) == 1:    
                    s_limits = np.array([s[:,0].min(), s[:,1].max(), s[:,0].max(), s[:,1].min()]).astype(int)
                     
                    if len(np.unique(np.append(stem_limits, s_limits))) == len(np.unique(s_limits)) + len(np.unique(stem_limits)):
                        if (abs(stem_limits[0]-s_limits[2]) == 1 or abs(stem_limits[1]-s_limits[3]) == 1 or
                            abs(s_limits[0]-stem_limits[2]) == 1 or abs(s_limits[1]-stem_limits[3]) == 1):
                            same_stem = np.append(same_stem, s, axis=0)
        same_stem = same_stem[np.argsort(same_stem[:,0] ),:]                    
        limits = np.array([same_stem[:,0].min(), same_stem[:,1].max(), same_stem[:,0].max(), same_stem[:,1].min()]).astype(int)
        n_stem = len(same_stem)
        param_ds = []
        t_stems = [tuple(s) for s in stem]

        for t_pair in t_stems:
            for pi in param_bp + param_stack + param_wc:
                if tuple(pi[1:3]) == t_pair or tuple(pi[1:3][::-1]) == t_pair:
                    pi[3] *= 1.2**(n_stem-1)
            p1 = t_pair[0]       
            p2 = t_pair[1]
            key, dim  = np.where(same_stem==t_pair)
            key = key[0]
            if abs(p1-p2) > 2:
                if key < len(same_stem)-1:
                    pul = same_stem[key+1][0]
                    pur = same_stem[key+1][1]
                    factor = 0.8**(abs(p1-pul)+abs(p2-pur)-2)
                    param_ang = np.append(param_ang, [[2, p1, p2, pul, pur, factor*n_stem * secon.k_ang, secon.angle]], axis=0)
        for pi in param_wc:
            if pi[1:3] in stem or pi[1:3][::-1] in stem:
                param_ds.append(pi)
        param_wc_ds.append(param_ds)

    sorted_params = np.empty((0,5))
    n_excl = 0
    for i in range(n):
        i_list = []
        for p in param_wc+param_bp+param_stack:
            if i == max(p[1:3]) and p[3] == 0:
                n_excl += 1
                continue
    #        elif i == max(p[1:3]) and abs(p[1]-p[2])<2:
    #            n_excl += 1
            elif i == max(p[1:3]):
                i_list.append(p)
        if len(i_list) == 0:
            continue
        diff = np.array([abs(p[1]-p[2]) for p in i_list])
        keys = np.argsort(diff)
        for k in keys:
            sorted_params = np.append(sorted_params, [i_list[k]], axis=0)
    assert (len(sorted_params)+n_excl == len(param_wc+param_bp+param_stack))        

    param_angle_180 = np.empty((0, 6))
    param_bulge = np.empty((0, 6))
    param_bulge_rep = np.empty((0, 6))
    paired = True
    pair_now = []
    for i in range(1, n-1):
        param_angle_180 = np.append(param_angle_180, [[8, i-1, i, i+1, secon.k_angle_straight, 0]], axis=0)        
    
    param_rep = np.empty((0, 5))
    param_rep_lr = np.empty((0, 5))
    i1, i2 = np.triu_indices(n, k=1)
    n_ij = i1.shape[0]
    _type = np.full((n_ij), 1)
    _k = np.full((n_ij), secon.k_rep2)
    _d = np.full((n_ij), secon.d_rep2)
    param_rep = np.append(param_rep, np.column_stack((_type, i1, i2, _k, _d)), axis=0)
    _type = np.full((n_ij), 7)
    _k = np.full((n_ij), secon.k_rep_lr)
    _d = np.full((n_ij), 0)
    param_rep_lr = np.append(param_rep_lr, np.column_stack((_type, i1, i2, _k, _d)), axis=0)

    return param_seq, param_rep, param_rep_lr, sorted_params, param_stem, param_ang, param_bulge, param_bulge_rep, param_angle_180   

def get_par(sorted_params, pos):
    i1 = sorted_params[:,1].astype(int)
    i2 = sorted_params[:,2].astype(int)
    d = np.linalg.norm(pos[i1]-pos[i2], axis=1)

    ii = np.argmin(d)
    return sorted_params[ii], ii
