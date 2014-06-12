from scipy.spatial import distance
import numpy as N
import pdbreader as pb
import tools as t


import itertools

def do_idx(j,l,b):
    # no bulges
    id = [N.arange(j,j+l).tolist()]

    if(b>0):
        # one bulge
        for i in xrange(j,j+l):
            t = N.arange(j,j+l+1).tolist()
            t.remove(i)
            id.append(t)
    if(b>1):
        #    two bulges
        for i in xrange(j,j+l+2):
            for k in xrange(i+1,j+l+2):
                t = N.arange(j,j+l+2).tolist()
                t.remove(i)
                t.remove(k)
                if t not in id:
                    id.append(t)
    return id
    
def ds_motif(args,files):

    fh = open(args.name,'w')

    print "# Finding 3D motifs..."
    pb.write_args(args,fh)
    l1 = args.l1
    l2 = args.l2

    atoms,ref_sequence = pb.get_coord(files[0])
    lcs,origo = t.coord2lcs(atoms[0])

    # calculate center of mass distances
    # this will be used to prune the search!
    delta_com = (N.sum(origo[0:l1],axis=0)/l1) - (N.sum(origo[l1:],axis=0)/l2)
    dd= N.sqrt(sum(delta_com**2))

    if(args.type=='modulus'):

        ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
        assert(ref_mat.shape[0] == args.l1+args.l2)

        ref_mat1 = ref_mat[0:l1,0:l1]
        ref_mat2 = ref_mat[l1:,l1:]
    

        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])
            
            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                flat_mat = mat.reshape(-1)
                
                # no bulges
                idx1 = []
                com1 = []
                for j in range(0,mat.shape[0]-l1+1):
                    tmp_idxs = do_idx(j,l1,args.bulges)
                    for tmp_idx in tmp_idxs:
                        if(tmp_idx[-1]>=mat.shape[0] or tmp_idx in idx1):
                            continue
                        red_mat1 = mat[tmp_idx,:][:,tmp_idx]
                        ermsd1 = t.calc_dist_1d(ref_mat1,red_mat1)
                        if(ermsd1 < args.treshold):
                            #for uu in tmp_idx:
                            #    print sequence[jj][uu],
                            #print ermsd1
                            idx1.append(tmp_idx)
                            com1.append(N.sum(origo[tmp_idx],axis=0)/l1)

                idx2 = []
                com2 = []
                for j in range(0,mat.shape[0]-l2+1):
                    tmp_idxs = do_idx(j,l2,args.bulges)
                    for tmp_idx in tmp_idxs:
                        if(tmp_idx[-1]>=mat.shape[0] or tmp_idx in idx2 ):
                            continue
                        red_mat2 = mat[tmp_idx,:][:,tmp_idx]
                        ermsd2 = t.calc_dist_1d(ref_mat2,red_mat2)
                        if(ermsd2 < args.treshold):
                            idx2.append(tmp_idx)
                            com2.append(N.sum(origo[tmp_idx],axis=0)/l2)

                        

                if(len(com1)==0 or len(com2) == 0 ):
                    continue

                dmat = distance.cdist(com1,com2)
                c_idx = (dmat<4.0*dd).nonzero()
                #print len(c_idx[0])
                for idx in range(len(c_idx[0])):
                    # skip overlapping
                    if(len(N.union1d(idx1[c_idx[0][idx]],idx2[c_idx[1][idx]]))!=l1+l2):
                        continue
                    all_idx = idx1[c_idx[0][idx]]+idx2[c_idx[1][idx]]
                    all_idx1 = idx1[c_idx[0][idx]]
                    all_idx2 = idx2[c_idx[1][idx]]
                    red_mat = mat[all_idx,:][:,all_idx]
                    ermsd = t.calc_dist_1d(ref_mat,red_mat)
                    if(ermsd < args.treshold):
                        seq = '; '
                        for el in all_idx1:
                            seq += (sequence[jj][el] + ' ')
                        seq += "* "
                        for el in all_idx2:
                            seq += (sequence[jj][el] + ' ')

                        string = '%8.5f %s %i %s \n' % (ermsd,files[ii],jj,seq)
                        fh.write(string)

                    

    if(args.type=='vector'):

        ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
        assert(ref_mat.shape[0] == args.l1+args.l2)

        ref_mat1 = ref_mat[0:l1,0:l1]
        ref_mat2 = ref_mat[l1:,l1:]
    

        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])
            
            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                
                # no bulges
                idx1 = []
                com1 = []
                for j in range(0,mat.shape[0]-l1+1):
                    tmp_idxs = do_idx(j,l1,args.bulges)
                    for tmp_idx in tmp_idxs:
                        if(tmp_idx[-1]>=mat.shape[0] or tmp_idx in idx1):
                            continue
                        red_mat1 = mat[tmp_idx,:][:,tmp_idx]
                        ermsd1 = t.calc_dist_nd(ref_mat1,red_mat1)
                        if(ermsd1 < args.treshold):
                            idx1.append(tmp_idx)
                            com1.append(N.sum(origo[tmp_idx],axis=0)/l1)

                idx2 = []
                com2 = []
                for j in range(0,mat.shape[0]-l2+1):
                    tmp_idxs = do_idx(j,l2,args.bulges)
                    for tmp_idx in tmp_idxs:
                        if(tmp_idx[-1]>=mat.shape[0] or tmp_idx in idx2 ):
                            continue
                        red_mat2 = mat[tmp_idx,:][:,tmp_idx]
                        ermsd2 = t.calc_dist_nd(ref_mat2,red_mat2)
                        if(ermsd2 < args.treshold):
                            idx2.append(tmp_idx)
                            com2.append(N.sum(origo[tmp_idx],axis=0)/l2)

                        

                if(len(com1)==0 or len(com2) == 0 ):
                    continue

                dmat = distance.cdist(com1,com2)
                c_idx = (dmat<4.0*dd).nonzero()
                print len(c_idx[0])
                for idx in range(len(c_idx[0])):
                    # skip overlapping
                    if(len(N.union1d(idx1[c_idx[0][idx]],idx2[c_idx[1][idx]]))!=l1+l2):
                        continue
                    all_idx = idx1[c_idx[0][idx]]+idx2[c_idx[1][idx]]
                    all_idx1 = idx1[c_idx[0][idx]]
                    all_idx2 = idx2[c_idx[1][idx]]

                    red_mat = mat[all_idx,:][:,all_idx]
                    ermsd = t.calc_dist_nd(ref_mat,red_mat)
                    if(ermsd < args.treshold):
                        seq = '; '
                        for el in all_idx1:
                            seq += (sequence[jj][el] + ' ')
                        seq += "* "
                        for el in all_idx2:
                            seq += (sequence[jj][el] + ' ')
                        string = '%8.5f %s %i %s \n' % (ermsd,files[ii],jj,seq)
                        fh.write(string)


####################### MOTIF ########################
