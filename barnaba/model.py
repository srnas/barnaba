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

import numpy as np
from scipy.spatial import distance
import sys
from . import residue as rd
from . import definitions


class Model:

    def __init__(self,data,set_lcs=False):

        ll = len(data)
        if(ll==0):
            sys.stderr.write("# Fatal error. No valid residues in pdb file \n")
            sys.exit(1)

        # array with all coordinates
        self.coords = np.array([data[i][5:8] for i in range(ll)])

        # lst with sequence
        self.sequence = []
        self.sequence_id = []
        self.chain = []
        self.residues = []
        resi_tmp = []

        for i in range(ll):
            if(i==ll-1):
                next_resid = "XXX"
            else:
                next_resid = data[i+1][9]
            resi_tmp.append(data[i])
            if(next_resid != data[i][9]):
                self.residues.append(rd.Residue(resi_tmp,i))
                self.sequence.append(data[i][8])
                self.sequence_id.append(data[i][9])
                resi_tmp = []

        
    def __getitem__(self, ind):
        return self.residues[ind]

    def set_coords(self,data):
        self.coords = data

    # Find index of backbone atoms 
    def set_bb_index(self):

        # handle exceptions
        self.missing = [False]*6*len(self.sequence)
        # first and last are always missing
        self.missing[0] = True
        self.missing[-1] = True
        self.missing[-2] = True
        self.chi_missing = [False]*len(self.sequence)
        
        bb_idx_tmp = []
        chi_idx = []
        for i in range(len(self.sequence)):
            residue = self.residues[i]
            atoms = [residue.get_idx(at_type) for at_type in definitions.rna_backbone]
            for jj in range(len(atoms)):
                # flag when meissing
                if(atoms[jj]==None):
                    nos = [x for x in range(6*i+jj-3,6*i+jj+1) if x>-1]
                    for el in nos:
                        self.missing[el+1] = True
            bb_idx_tmp.extend(atoms)
            
            # now for chi angle
            atoms = [residue.get_idx("O4'"),residue.get_idx("C1'")]
            if(self.sequence[i] in definitions.pur):
                atoms.append(residue.get_idx("N9"))
                atoms.append(residue.get_idx("C4"))
            if(self.sequence[i] in definitions.pyr):
                atoms.append(residue.get_idx("N1"))
                atoms.append(residue.get_idx("C2"))
            if(None not in atoms):
                chi_idx.append(atoms)
            else:
                self.missing[i] = True
                

        # handle missing atoms
        bb_idx = []
        for ii in range(len(bb_idx_tmp)-3):
            idxs = bb_idx_tmp[ii:ii+4]
            if(None not in idxs):
                bb_idx.append(idxs)
        bb_angles = self.bb_idx = np.array(bb_idx)
        chi_angles = self.chi_idx = np.array(chi_idx)

        



    def get_mat_score(self,cutoff,ii=[]):

        ll = len(ii)
        if(ll==0):
            ll = len(self.sequence)
            ii = range(ll)            
        
        dotp,m_idx = self.get_mat_score(cutoff,ii)
        dotp_scale = dotp*np.array(definitions.scale)[np.newaxis,:]

        # set to zero when norm is larger than cutoff
        dotp_norm = np.sqrt(np.sum(dotp_scale**2,axis=1))
        dotp[dotp_norm>cutoff] = 0.0
        nonzero = dotp[~np.all(dotp == 0, axis=1)]
        
        return nonzero.T

    

    def string_pdb(self,idx=None,noP=False,center=False):

        if(idx==None):
            idx = range(len(self.residues))

        string = "MODEL \n"
        string += "REMARK "
        string += " ".join([self.sequence_id[i] for i in idx])
        string += "\n"

        com = np.zeros(3)
        if(center):
            # calculate center of mass
            com =  np.average(np.array([np.average(self.coords[self.residues[i].first:self.residues[i].last+1],axis=0) for i in idx]),axis=0)

        for serial,i in enumerate(idx):
            coords = self.coords[self.residues[i].first:self.residues[i].last+1]-com
            if(noP==True and serial==0):
                string += self.residues[i].pdb_string(coords,noP=True)
            else:
                string += self.residues[i].pdb_string(coords)

        string += "ENDMDL\n"
        return string





