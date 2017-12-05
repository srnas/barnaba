import mdtraj as md
import functions
import nucleic
import numpy as np
import kde as kde
import sys

class Escore:

    def __init__(self,references,cutoff=1.58,bandwidth=0.25):

        mats = []
        for i in range(len(references)):
            pdb = md.load(references[i])    
            nn_ref = nucleic.Nucleic(pdb.topology,modified=False)
            coords = pdb.xyz[0,nn_ref.indeces_lcs]
            mat = functions.calc_scoremat(coords,cutoff)
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

        print "AAA"
        print sample, topology
        if(topology==None):
            traj = md.load(sample)
        else:
            traj = md.load(sample,top=topology)
        warn = "# Loaded sample %s \n" % sample
        sys.stderr.write(warn)
        
        nn = nucleic.Nucleic(traj.topology,modified=False)
        scores = []
        for j in xrange(traj.n_frames):
            coords = traj.xyz[j,nn.indeces_lcs]
            mat = functions.calc_scoremat(coords,self.cutoff+0.2)
            scores.append(np.sum(self.kernel(10.0*mat)))
        return scores
