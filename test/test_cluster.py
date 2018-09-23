from __future__ import absolute_import, division, print_function
import sys
import numpy as np
import os
from comp_mine import comp

import barnaba as bb
import barnaba.cluster as cc
import mdtraj as md

#import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set_style("white")

cwd = os.getcwd()
outdir = "%s/test/tmp" % cwd
refdir = "%s/test/reference/" % cwd
os.system("mkdir -p %s" % (outdir))

top = "%s/test/data/UUCG.pdb" % cwd
traj = "%s/test/data/UUCG.xtc" % cwd

def test_cluster():
    # first, calculate all g-vectors
    print("# Calculating G-vectors")
    gvec,seq = bb.dump_gvec(traj,top)
    lent = gvec.shape[0]
    gvec = gvec.reshape(lent,-1)[::5]
    
    print("# Calculating PCA. gvec shape: ", gvec.shape)
    # calculate PCA
    v,w = cc.pca(gvec,nevecs=3)
    print("# Cumulative explained variance of component: 1=%5.1f 2:=%5.1f 3=%5.1f" % (v[0]*100,v[1]*100,v[2]*100))
    print("# DBSCAN clustering...")
    # do DBSCAN clustering. eps and min_samples need to be adjusted.
    new_labels, center_idx = cc.dbscan(gvec,range(gvec.shape[0]),eps=0.6/np.sqrt(8.),min_samples=10)
    print("DONE!")

    # create color palette. gray, small points for unassigned clusters.
    #cp = sns.color_palette("hls",len(center_idx)+1)
    #colors = [cp[j-1] if(j!=0) else (0.77,0.77,0.77) for j in new_labels]
    size = [5 if(j!=0) else 0.25 for j in new_labels]
    # do scatterplot
    #plt.scatter(w[:,0],w[:,1],s=size,c=colors)

    # now dump centroids and print labels on plot
    
    print("# Dump PDB centroids")
    t = md.load(traj, top=top)
    idxs = [ii for ii,kk in enumerate(new_labels) if(kk==0)]
    for i,k in enumerate(center_idx):
        t[k].save_pdb("%s/cluster_%03d.test.pdb" % (outdir,i))
        comp("%s/cluster_%03d.test.pdb" % (refdir,i))
        #plt.text(w[k,0],w[k,1],str(i),ha='center',va='center')
        idxs = [ii for ii,kk in enumerate(new_labels) if(kk==i+1)]

    # save figure
    #plt.savefig("clusters.png",dpi=300)
    #plt.close()


