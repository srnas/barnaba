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

""" scoring function for RNA structure prediction """

from __future__ import absolute_import, division, print_function

# Make sure that range returns an iterator also in python2 (using future module)
from builtins import range

import mdtraj as md
import numpy as np
import sys
from . import kde
from . import calc_mats as ff
from . import nucleic

class Escore:

    """ Constructor """
    def __init__(self,references,cutoff=1.58,bandwidth=0.25):

        mats = []
        for i in range(len(references)):
            pdb = md.load(references[i])    
            nn_ref = nucleic.Nucleic(pdb.topology,modified=False)
            coords = pdb.xyz[0,nn_ref.indeces_lcs]
            mat = ff.calc_scoremat(coords,cutoff)
            mats.extend(mat)
        mats = np.asarray(mats)
        # set kernel density
        self.kernel = kde.gaussian_kde(10.0*mats)
        self.cutoff = cutoff
        # set kernel to 0.25(def)
        self.kernel.set_bandwidth(bandwidth)
        warn =  "# KDE computed. Bandwidth= %5.2f " % (self.kernel.factor)
        warn += " using %d base-pairs" % mats.shape[1]
        sys.stderr.write(warn)
        
    def score(self,sample,topology=None):
        
        """ Score """
        if(topology==None):
            traj = md.load(sample)
        else:
            traj = md.load(sample,top=topology)
        warn = "# Loaded sample %s \n" % sample
        sys.stderr.write(warn)
        
        nn = nucleic.Nucleic(traj.topology,modified=False)
        scores = []
        for j in range(traj.n_frames):
            coords = traj.xyz[j,nn.indeces_lcs]
            mat = ff.calc_scoremat(coords,self.cutoff+0.2)
            scores.append(np.sum(self.kernel(10.0*mat)))
        return scores
