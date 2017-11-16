#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014 Sandro Bottaro (sbottaro@sissa.it) and Giovanni Pinamonti (giopina@sissa.it)

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
from scipy.linalg import eigh
from scipy.spatial import distance
import mdtraj as md
import definitions
import nucleic

class Annotate:

    def __init__(self,pdb):

        # load pdb 
        cur_pdb = md.load_pdb(pdb)
        topology = cur_pdb.topology

        nn = nucleic.Nucleic(topology)
        for ii in nn.ok_residues:
            print ii
