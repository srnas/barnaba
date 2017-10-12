import sys
import nucleic
import numpy as np
import tools as tools
import mdtraj as md

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




