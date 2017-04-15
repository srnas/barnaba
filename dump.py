#import sys
import nucleic
import numpy as np

def get_string(seq,data,hread=True):

    s = ""
    if(hread==False):
        data = data.reshape(-1)
        #s = "%10s " % (fname)
        s += "".join([("%15e " % el) for el in data])
        return s + "\n"
    else:
        for j in range(data.shape[0]):
            for k in range(data.shape[1]):
                if(np.sum((data[j,k])**2) > 0.0):
                    s += '%15s %15s ' % (seq[j],seq[k])
                    s += "".join(['%15e ' % data[j,k,l] for l in range(len(data[j,k]))])
                    s+= "\n"
        return s



def dump_rvec(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    rvecs = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        rvecs.append(nn.get_rmat(coords_lcs,cutoff))

    return nn.rna_seq,np.asarray(rvecs)


def dump_gvec(traj,cutoff=2.4):
        
    top = traj.topology
    nn = nucleic.Nucleic(top)
    gvecs = []
    for i in xrange(traj.n_frames):
        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        gvecs.append(nn.get_gmat(coords_lcs,cutoff))
    return nn.rna_seq,np.asarray(gvecs)
