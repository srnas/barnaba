#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014 Sandro Bottaro (sbottaro@sissa.it)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#import reader as reader
import numpy as np
import mdtraj as md
import btools as bt
from scipy.spatial.distance import cdist,squareform
from scipy.linalg import eig, eigh
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics.pairwise import pairwise_distances

def cluster(args):
    
    print "# Calculating ..."
    fh = open(args.gvec)

    gvecs = []
    labels = []
    for line in fh:
        if("#" not in line):
            gvecs.append([float(x) for x in line.split()[1:]])
            labels.append(line.split()[0])
    fh.close()
    gvecs = np.array(gvecs)
    nn = 6
    print "# read ", gvecs.shape, "gmats"
    slen = np.sqrt(gvecs.shape[1]/4)
    if((slen).is_integer()):
        print "# Sequence lenght =",slen
    else:
        print "# Error - check your gmat file"
        sys.exit(1)

    #fh = open(args.name,'a')
    
    # kernelPCA
    nc = args.ncluster
    dmat= None
    ll = gvecs.shape[0]

    if(args.alg=="kernel"):

        print "# Calculate pairwise distance"
        dmat = pairwise_distances(gvecs,n_jobs=-1)/np.sqrt(float(slen))
        # adjacency matrix
        print "# done"
        kmat= np.exp(-0.5*(dmat*dmat)/(args.eps*args.eps))

        # set diagonal to zero
        kmat[np.arange(ll),np.arange(ll)] = 0.0
        # inverse degree matrix
        degree_sinv = np.sqrt(1./np.sum(kmat,axis=1))
        kmat_sym = np.outer(degree_sinv,degree_sinv)*kmat

        #kmat_sym = np.identity(dmat.shape[0]) - np.outer(degree_sinv,degree_sinv)*kmat
        
        #kmat_norm = np.array([kmat[i,:]/degree_s for i in range(kmat.shape[0])])
        print "# diagonalize"
        v,w = np.linalg.eigh(kmat_sym)
        idx2 = (np.abs(v)).argsort()[::-1]
        v = v[idx2]
        w = w[:,idx2]
        kmeans_w = np.copy(w[:,1:nc])

        # normalize
        for ii in range(kmeans_w.shape[0]):
            kmeans_w[ii] /= np.sqrt(np.sum(kmeans_w[ii]**2))
            
        # do kmeans
        k_means = KMeans(init='k-means++',n_init=50,n_clusters=nc,random_state=144)
        k_means.fit(kmeans_w)
        
        cluster_labels = k_means.labels_
        centers = k_means.cluster_centers_
    
        # create array where each row contains members
        cluster_members = [[] for x in range(nc)]
        for j in range(kmeans_w.shape[0]):
            cluster_members[cluster_labels[j]].append(j)

    if(args.alg=="DBSCAN"):
        ns = int(gvecs.shape[0]/200)
        #ns = 50
        print "# Running dbscan: eps=%f, min_samples=%d" % (args.eps,ns)
        db = DBSCAN(eps=args.eps, min_samples=50).fit(gvecs)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        cluster_labels = db.labels_
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
        cluster_members = [[] for x in range(n_clusters_+1)]
        for j in range(gvecs.shape[0]):
            #if(cluster_labels[j] != -1):
            cluster_members[cluster_labels[j]].append(j)
        
        nc = n_clusters_+1


        
    if(args.mode=="PCA"):

        data =np.copy(gvecs)
        
        # subtract average
        avg  = []
        for j in range(data.shape[1]):
            data[:,j] -= np.average(data[:,j])

        # calculate covariance
        cov = np.cov(data.T)
        # diagonalize
        vv, ww = np.linalg.eigh(cov)
        idx2 = (np.abs(vv)).argsort()[::-1]
        evals = vv[idx2]
        evecs = ww[:,idx2]
        
        sum_evals = np.sum(evals)
        v = [np.sum(evals[:jj])/sum_evals for jj in range(1,len(evals)+1)]
            
        # pad with ones...
        v += [1.0]*(data.shape[0]-len(v))
        w = []
        for c in range(nn):
            w.append(np.dot(evecs[:,c],data.T))
        w = np.array(w).T

    if(args.mode=="kernel"):

        print "# Calculate pairwise distance"
        ll = gvecs.shape[0]
        if(dmat==None): dmat = pairwise_distances(gvecs,n_jobs=-1)/np.sqrt(float(slen)) 

        # adjacency matrix
        print "# done"
        kmat= np.exp(-0.5*(dmat*dmat)/(args.eps_p*args.eps_p))
        # set diagonal to zero
        kmat[np.arange(ll),np.arange(ll)] = 0.0
        # inverse degree matrix
        degree_sinv = np.sqrt(1./np.sum(kmat,axis=1))
        kmat_sym = np.outer(degree_sinv,degree_sinv)*kmat

        #kmat_sym = np.identity(dmat.shape[0]) - np.outer(degree_sinv,degree_sinv)*kmat
        
        #kmat_norm = np.array([kmat[i,:]/degree_s for i in range(kmat.shape[0])])
        print "# diagonalize"
        v,w = np.linalg.eigh(kmat_sym)
        idx2 = (np.abs(v)).argsort()[::-1]
        v = v[idx2]
        w = w[:,idx2]
        
        
        

    fh_c = open(args.name + ".centers.dat",'w')
    ss = "# %2s %70s %4s %7s %7s %7s \n" % ("n","label","size","error","max d", "avg dist")
    for j in range(nc):
                
        # get all inter-cluster distances
        gvec_tmp = np.array([gvecs[k] for k in cluster_members[j]])
        labels_tmp = [labels[k] for k in cluster_members[j]]               
        dist = cdist(gvec_tmp,gvec_tmp)/np.sqrt(slen)
        
        # find the one that minimises the sum from all others
        dd = np.sum(dist,axis=1)
        min_idx = np.argmin(dd)
        idx_centroid = labels.index(labels_tmp[min_idx])
        # calculate distance between centroid and geometrical center 
        try:
            err = kmeans_w[cluster_members[j][min_idx]]-np.array(centers[j])                
            err = np.sqrt(np.sum(err**2))
        except:
            err = 0.0
        ss += "%4d %70s %5d %7.3f %7.3f %7.3f \n" % \
            (j,labels_tmp[min_idx], len(cluster_members[j]),err, np.max(dist),np.average(dist) )
    fh_c.write(ss)
    fh_c.close()

            
    fh = open(args.name + ".cluster.dat",'w')
    header = "#%69s %14s " % ("Label","Eigenvalue")

    for j in range(nn):
        header += "%14s %5s \n" % ("Eigenvector","Cluster")
    for j in range(len(v)):
        # eigenvalues
        stri = "%70s %14.5e " % (labels[j],float(v[j]))
        for k in range(nn):
            stri += "%14.5e " % (float(w[j,k]))
        stri += "%5d \n" % (cluster_labels[j])
        fh.write(stri)
    fh.close()  


    return 0


