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
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from scipy.spatial.distance import squareform,pdist

def pca(gvecs,nevecs=6,sample_weight=None):

    # subtract average
    data = gvecs-np.average(gvecs,axis=0,weights=sample_weight)[np.newaxis,:]       

    # calculate covariance
    cov = np.cov(data.T,aweights=sample_weight)
    # diagonalize
    vv, ww = np.linalg.eigh(cov)
    idx2 = (np.abs(vv)).argsort()[::-1]
    evals = vv[idx2]
    evecs = ww[:,idx2]
    sum_evals = np.sum(evals)
    v = [np.sum(evals[:jj])/sum_evals for jj in range(1,len(evals)+1)]
    
    w = []
    for c in range(nevecs):
        w.append(np.dot(evecs[:,c],data.T))
    w = np.array(w).T
    return v,w


def dbscan(gvecs,labels,eps,min_samples,sample_weight=None):

    
    slen = np.sqrt(gvecs.shape[1]/4.)
    #print(slen)
    eps *=np.sqrt(slen)
    db = DBSCAN(eps=eps, min_samples=min_samples,algorithm='brute').fit(gvecs,sample_weight=sample_weight)
    #db = DBSCAN(eps=eps, min_samples=min_samples).fit(gvecs)
    cluster_labels = db.labels_
    
    n_clusters_ = len(set(cluster_labels))
    print('# eps:%5.3f min_samples:%d  nclusters: %d' % (eps,min_samples, n_clusters_-1))
    print("#  silhouette score: %0.4f"
          % metrics.silhouette_score(gvecs, cluster_labels)),
    non_l = np.where(cluster_labels!=-1)
    print("# Avg silhouette: %0.4f " % (np.average(metrics.silhouette_samples(gvecs, cluster_labels)[non_l])))
    print("# assigned samples :%d total samples:%d " % (len(non_l[0]), len(cluster_labels)))
        
    cluster_members = [[] for x in range(n_clusters_-1)]

    for j in range(gvecs.shape[0]):
        ii1 = cluster_labels[j]
        if(ii1!=-1):
            cluster_members[ii1].append(j)

    scenters = np.argsort([len(x) for x in cluster_members])[::-1]
    print("# %2s %4s %20s %20s %20s %20s %s " % ("N","size","max eRMSD (IC)","med eRMSD (IC)","max eRMSD (centroid)","med eRMSD (centroid)","center"))

    labels_map = [None]*n_clusters_
    labels_map[-1] = 0
    center_idx = []
    for o,ii1 in enumerate(scenters):
        
        gvec_tmp = gvecs[cluster_members[ii1]]
        dists = squareform(pdist(gvec_tmp))/np.sqrt(float(slen))
        labels_tmp = [labels[k] for k in cluster_members[ii1]]               
        dd = np.sum(dists,axis=1)
        min_idx = np.argmin(dd)
        print("# %02d %04d %20.3f %20.3f %20.3f %20.3f %02d %s" % (o,len(cluster_members[ii1]), np.max(dists), np.median(dists),np.max(dists[min_idx]),np.median(dists[min_idx]),ii1,labels_tmp[min_idx]))
        labels_map[ii1] = o+1
        center_idx.append(labels.index(labels_tmp[min_idx]))
    new_labels = [labels_map[cluster_labels[k]] for k in range(len(cluster_labels))]
    return new_labels, center_idx
        
