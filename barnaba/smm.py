from __future__ import absolute_import, division, print_function
import sys
import numpy as np
from scipy.spatial.distance import cdist,squareform
import argparse
from sklearn.cluster import KMeans
import bisect

class SMM:


    def __init__(self,gvecs,eps,weights=[]):

        if(len(weights)==0):
            weights = np.ones(gvecs.shape[0])
        else:
            print("#not fully tested w non-uniform weights")
            weights = np.array(weights)
            assert weights.shape[0] == gvecs.shape[0]
            
        # infer lenght        
        slen = np.sqrt(gvecs.shape[1]/4)
        if((slen).is_integer()):
            print("# Sequence lenght =",slen)
            self.slen = slen
        else:
            print("# Error - check your gmat file")
            sys.exit(1)

        # calculate distance matrix
        dmat_sq = cdist(gvecs,gvecs,'sqeuclidean')
        small = (dmat_sq<0.001).nonzero()
        print("# Number of identical elements...",len(small[0])/2-self.slen)

        kmat= np.exp(-(0.5*dmat_sq)/(eps*eps*self.slen))
        
        print("# Iterative normalization")
        tol = 1000
        tolerance = 1.0E-10
        it = 1
        rang = range(kmat.shape[0])
        while tol > tolerance:

            degree_s = np.asarray([np.sum(weights*kmat[i,:]) for i in rang])
            diff1 = np.abs(1.0 - np.max(degree_s))
            diff2 = np.abs(1.0 - np.min(degree_s))
            tol = max(diff1,diff2)
            degree_sinv = 1./np.sqrt(degree_s)
            sqmat = np.outer(degree_sinv,degree_sinv) 
            kmat *= sqmat
            it += 1
        print("# Done", it , "iterations")

        mat = np.multiply(kmat,weights[:,np.newaxis])

        print("# calculated transition matrix mat")
        



