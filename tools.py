import definitions
import itertools as its
import numpy as np
import re
import sys

def get_pattern(query):
    # build pattern for regular expression
    pattern = ""
    for res in query:
        assert res in definitions.known_abbrev, "# Fatal error: character %s not known. Use AUCG/NYR" % (res)
        if(res in definitions.rna):
            pattern += res
        else:
            if(res == "N"):
                pattern += "[AUCG]"
            if(res == "Y"):
                pattern += "[UC]"
            if(res == "R"):
                pattern += "[AG]"
    return pattern



def get_idx(sequence,query,bulges=0):
    
    ll = len(query)
    seq_str  = "".join(sequence)
    pattern = get_pattern(query)

    # generate first all possible indeces
    indeces = []
    for b in range(bulges+1):

        # create position of insertions
        comb= its.combinations(range(1,ll+b-1),b)
        for it1 in comb:
            idx1 = range(ll+b)
            # remove them from list
            for it2 in it1:
                idx1.remove(it2)
            for start in range(len(sequence)-ll-b+1):
                # get sequence
                idx2 =  [idx1[i] + start for i in range(ll)]
                substr = "".join([sequence[i] for i in idx2])                
                if(re.match(pattern,substr) != None):
                    indeces.append(idx2)
    return indeces

def is_complementary(r1,r2):

    if(r1=="A" and r2 == "U"): return True
    if(r2=="A" and r1 == "U"): return True
    if(r1=="G" and r2 == "C"): return True
    if(r2=="G" and r1 == "C"): return True
    return False

def wc_gaussian(vec):

     dev1  = vec - definitions.wc_mean
     maha1 = np.einsum('...k,...kl,...l->...', dev1, definitions.inv_sigma, dev1)
     return (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha1)
     
            
        
def chain_consistency(indeces,seq_id):

    removeme = []
    for ii in xrange(len(indeces)):
        idx = indeces[ii]
        chains = [x.split("_")[-1] for x in seq_id[idx[0]:idx[-1]+1]]
        if (any(x!=chains[0] for x in chains)):
            removeme.append(ii)
            continue
        nums = [int(x.split("_")[0]) for x in seq_id[idx[0]:idx[-1]+1]]
        for i2 in xrange(1,len(nums)):
            if(nums[i2]!=(nums[i2-1]+1)):
                removeme.append(ii)
                continue
    # remove 
    for el in removeme[::-1]:
        indeces.pop(el)
    #return indeces


def dihedral(vecs):

    # difference vectors b0 is reversed
    b0 = vecs[:,1] - vecs[:,0]
    b1 = vecs[:,1] - vecs[:,2]
    b2 = vecs[:,2] - vecs[:,3]
    # norm
    norm_sq = np.sum(b1**2,axis=1)
    norm_sq_inv = 1.0/norm_sq
    large_norm =  (norm_sq>2.9).nonzero()
    if(len(large_norm[0])>0):
        err = "# Warning: a bond lenght is suspiciously large % \n"
        print norm_sq[norm_sq>2.9]
        print large_norm
        exit()
        norm_sq_inv[large_norm] = float('nan')
        sys.stderr.write(err)

    #print (np.sum(b0*b1,axis=1)*b1).shape
    v0 = b0 - b1*((np.sum(b0*b1,axis=1)*norm_sq_inv)[:,np.newaxis])
    v2 = b2 - b1*((np.sum(b0*b2,axis=1)*norm_sq_inv)[:,np.newaxis])
    x = np.sum(v0*v2,axis=1)
    m = np.cross(v0,b1)*np.sqrt(norm_sq_inv[:,np.newaxis])
    y = np.sum(m*v2,axis=1)
    return np.arctan2( y, x )

    

