import sys
import numpy as np
from scipy.spatial.distance import cdist,squareform
import argparse
from sklearn.cluster import KMeans
import bisect



def print_kmat(labels,v,w,nevec,filename,plus=[]):

    # Print to file
    fh = open(filename,'w')
    for j in range(len(v)):
        # eigenvalues
        stri = "%14.5e " % (float(v[j]))
        for k in range(nevec):
            stri += "%14.5e " % (float(w[j,k]))
        if(len(plus)==len(v)):
            stri += "%3d " % (plus[j])
        stri += "%s \n" % labels[j]
        fh.write(stri)
    fh.close()  

        

class SMM:
    
    def gvec2dmat(self):

        self.dmat = cdist(self.gvecs,self.gvecs)/np.sqrt(self.slen)
        ww = (self.dmat<0.001).nonzero()
        print "# Number of identical elements...",len(ww[0])/2

    def calc_kmat(self,eps):

        dist_sq  = self.dmat*self.dmat
        kmat= np.exp(-0.5*(dist_sq)/(eps*eps))
        
        print "# Iterative normalization"
        tol = 1000
        tolerance = 1.0E-10
        it = 1
        while tol > tolerance:
            degree_s = np.sum(kmat,axis=1)
            diff1 = np.abs(1.0 - np.max(degree_s))
            diff2 = np.abs(1.0 - np.min(degree_s))
            tol = max(diff1,diff2)
            degree_sinv = 1./np.sqrt(degree_s)
            sqmat = np.outer(degree_sinv,degree_sinv) 
            kmat *= sqmat
            it += 1
        print "# Done", it , "iterations"
        return kmat

    

    def __init__(self,gvec_file,ref_file,eps,name):

        # initalize variables
        data = []
        self.labels = []
        self.name = name
        self.evec = None
        self.eval = None
        self.cluster_members = None
        
        # read reference file with source and sink
        fh = open(ref_file)
        for line in fh:
            if("#" not in line):
                data.append([float(x) for x in line.split()[1:]])
                ll = line.split()[0]
                self.labels.append(ll)
                if(len(data)==3):
                    print "# Fatal. Only two GVEC in reference file."
                    sys.exit()
        fh.close()
        self.idx_source = 0
        self.idx_sink = 1

        # read GVECTOR
        fh = open(gvec_file)
        for line in fh:
            if("#" not in line):
                vv = [float(x) for x in line.split()[1:]]
                if(len(vv)!=len(data[0])):
                    print "# skipping", line.split()[0]
                    continue
                data.append(vv)
                ll = line.split()[0]
                self.labels.append(ll)
        fh.close()
        data = np.array(data)
        self.gvecs = data
        print "# Read ", data.shape[0], "g-vectors (",data.shape,")"

        # infer lenght        
        slen = np.sqrt(data.shape[1]/4)
        if((slen).is_integer()):
            print "# Sequence lenght =",slen
            self.slen = slen
        else:
            print "# Error - check your gmat file"
            sys.exit(1)

        # calculate distance matrix
        self.gvec2dmat()
        # calculate transition matrix
        self.kmat = self.calc_kmat(eps)

            
    def do_trajectory(self,ntraj,distance):

        # start from A-form 
        cur_idx = self.idx_source
        
        inA = self.dmat[self.idx_source]<distance
        inB = self.dmat[self.idx_sink]<distance
        print "# Number of samples in A-form:", (len(inA.nonzero()[0]))
        print "# Number of samples in Tloop-form:", (len(inB.nonzero()[0]))
        
        # calculate cumulative
        cumsum = np.cumsum(self.kmat,axis=1)


        traj = []
        trajs = []
       
        fh_t = open(self.name + ".traj.dat",'w')
        
        while(len(trajs) < ntraj):
            cur_idx = bisect.bisect_right(cumsum[cur_idx], np.random.random())
            traj.append(cur_idx)
            if(len(traj)==1E+06): 
                print "# Length of the stochastic trajectory exceeding 1M steps. Perhaps the Network is not connected? "
                print "# .. ABORTING .."
                sys.exit(1)
            if(inA[cur_idx]):
                if(inB[traj[0]]):
                    trajs.append(traj[::-1])
                    fh_t.write(" ".join([str(x) for x in traj[::-1]]) + "\n")
                traj = [cur_idx]
            if(inB[cur_idx]):
                if(inA[traj[0]]):
                    trajs.append(traj)
                    fh_t.write(" ".join([str(x) for x in traj]) + "\n")
                traj = [cur_idx]

        fh_t.close()

    def calc_evecs(self):

        if(self.eval==None):
            v,w = np.linalg.eigh(self.kmat)
            idx2 = (np.abs(v)).argsort()[::-1]
            self.eval = v[idx2]
            self.evec = w[:,idx2]

        
    def calc_kmeans(self,nc):

        if(self.cluster_members ==None):

            # calculate eigenvalues and eigenvectors
            self.calc_evecs()
            
            # take first nc coordinates
            kmeans_w = np.copy(self.evec[:,1:nc])
            # normalize
            for ii in range(kmeans_w.shape[0]):
                kmeans_w[ii] /= np.sqrt(np.sum(kmeans_w[ii]**2))
            
            # do kmeans
            k_means = KMeans(init='k-means++',n_init=50,n_clusters=nc,random_state=144)
            k_means.fit(kmeans_w)
            
            self.cluster_labels = k_means.labels_
            centers = k_means.cluster_centers_

            # create array where each row contains members
            cluster_members = [[] for x in range(nc)]
            for j in range(kmeans_w.shape[0]):
                cluster_members[self.cluster_labels[j]].append(j)

            self.cluster_members = cluster_members

            fh_c = open(self.name + ".centers.dat",'w')
            ss = "# %2s %70s %4s %7s %7s %7s \n" % ("n","label","size","error","max d", "avg dist")
            for j in range(nc):
                
                # get all inter-cluster distances
                gvec_tmp = np.array([self.gvecs[k] for k in self.cluster_members[j]])
                labels_tmp = [self.labels[k] for k in self.cluster_members[j]]               
                dist = cdist(gvec_tmp,gvec_tmp)/np.sqrt(self.slen)
        
                # find the one that minimises the sum from all others
                dd = np.sum(dist,axis=1)
                min_idx = np.argmin(dd)
                idx_centroid = self.labels.index(labels_tmp[min_idx])

                # calculate distance between centroid and geometrical center 
                err = kmeans_w[self.cluster_members[j][min_idx]]-np.array(centers[j])                
                ss += "%4d %70s %5d %7.3f %7.3f %7.3f \n" % \
                    (j,labels_tmp[min_idx], len(self.cluster_members[j]),np.sqrt(sum(err**2)), np.max(dist),np.average(dist) )
            fh_c.write(ss)
            fh_c.close()

        
    def calc_flux(self):

        def get_map(cg,member):

            sort1 = [sorted(x) for x in cg]
            sort2 = [sorted(x) for x in member]

            mapping = [-1]*len(cg)
            ll = len(cg)
            for ii in range(ll):
                for jj in range(ll):
                    if(len(sort1[ii]) != len(sort2[jj])): continue

                    for el1,el2 in zip(sort1[ii],sort2[jj]):
                        if(el1 != el2):
                            break
                        mapping[ii] = jj
            return mapping


        # find cluster of helix and loop
        cluster_h = self.cluster_labels[self.idx_source]
        cluster_l = self.cluster_labels[self.idx_sink]
        if(cluster_h == cluster_l):
            print "# Fatal error. Source and sink belong to the same cluster (eRMSD= %f)" % (self.dmat[self.idx_source,self.idx_sink])
            sys.exit(1)
            
        A = self.cluster_members[cluster_h]
        B = self.cluster_members[cluster_l]
        
        # get fluxes
        # import pyemma
        import pyemma.msm as msm      

        M = msm.MSM(self.kmat)
        fluxAB = msm.tpt(M,A,B)
        cg, cgflux =  fluxAB.coarse_grain(self.cluster_members)
        paths, path_fluxes = cgflux.pathways(fraction=0.99)


        # find mapping (pyemma messes up the indeces)
        mapping = get_map(cg,self.cluster_members)

        # write to file
        fh_c = open(self.name + ".pathways.dat",'w')
        fh_c.write("# cluster source " + str(cluster_h) + "\n")
        fh_c.write("# cluster sink  " + str(cluster_l) + "\n")

        stri = "# Flux % - List of clusters \n"
        for i in range(len(paths)):
            stri += "%5.3f - " %  (path_fluxes[i] / np.sum(path_fluxes))
            for el in paths[i]:
                stri += "%3i " % mapping[el]
            stri += "\n"
        fh_c.write(stri)
        fh_c.close()

    def print_projections(self,mode,large_eps):
        

        # print first 5/6 coordinates only
        nn = 6
        
        if(mode=="kernel"):
            # calculate kmat with large eps
            kmat_large = self.calc_kmat(large_eps)
            # diagonalize
            v,w = np.linalg.eigh(kmat_large)
            idx2 = (np.abs(v)).argsort()[::-1]
            v = v[idx2]
            w = w[:,idx2]

        if(mode=="PCA"):
            data =np.copy(self.gvecs)
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
            
            
        fh = open(self.name + ".smm.dat",'w')
        header = "#%69s %14s " % ("Label","Eigenvalue")
        for j in range(nn):
            header += "%14s " % ("Eigenvector")
        if(self.cluster_members != None):
            header += "%5s" % ("Clust")
        fh.write(header + "\n")
        for j in range(len(v)):
            # eigenvalues
            stri = "%70s %14.5e " % (self.labels[j],float(v[j]))
            for k in range(nn):
                stri += "%14.5e " % (float(w[j,k]))
            # add cluster 
            if(self.cluster_members != None):
                stri += "%5d " % (self.cluster_labels[j])
            stri += "\n"
            fh.write(stri)
        fh.close()  

        

def smm(args):
    
    smm = SMM(args.gvec,args.ref,args.eps,args.name)
    
    smm.calc_kmeans(args.ncluster)
    smm.calc_flux()

    if(args.ntraj>0):
        smm.do_trajectory(args.ntraj,args.dist)




    # print labels. use larger eps for calculating projection 
    smm.print_projections(args.mode,args.epsc)



