from scipy.spatial import distance
import numpy as N
import pdbreader as pb
import tools as t


import itertools

def ds_motif(args,files):

    fh = open(args.name,'w')

    print "# Finding 3D motifs..."
    pb.write_args(args,fh)
    l1 = args.l1
    l2 = args.l2

    atoms,ref_sequence = pb.get_coord(files[0])
    lcs,origo = t.coord2lcs(atoms[0])
    ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
    assert(ref_mat.shape[0] == args.l1+args.l2)

    ref_mat1 = ref_mat[0:l1,0:l1]
    ref_mat2 = ref_mat[l1:,l1:]
    
    # calculate center of mass distances
    # this will be used to prune the search!
    delta_com = (N.sum(origo[0:l1],axis=0)/l1) - (N.sum(origo[l1:],axis=0)/l2)
    dd= N.sqrt(sum(delta_com**2))

    if(args.type=='modulus'):

        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])
            
            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                flat_mat = mat.reshape(-1)
                
                
                idx1 = []
                com1 = []
                for j in range(0,mat.shape[0]-l1+1):
                    tmp_idx = N.arange(j,j+l1).tolist()
                    red_mat1 = mat[tmp_idx,tmp_idx]
                    ermsd1 = t.calc_dist_1d(ref_mat1,red_mat1)
                    if(ermsd1 < args.treshold):
                        idx1.append(tmp_idx)
                        com1.append(N.sum(origo[tmp_idx],axis=0)/l1)

                idx2 = []
                com2 = []
                for j in range(0,mat.shape[0]-l2+1):
                    red_mat2 = mat[j:j+l2,j:j+l2]
                    ermsd2 = t.calc_dist_1d(ref_mat2,red_mat2)
                    if(ermsd2 < args.treshold):
                        idx2.append(N.arange(j,j+l2).tolist())
                        com2.append(N.sum(origo[j:j+l2],axis=0)/l2)
                        
                if(len(com1)==0 or len(com2) == 0 ):
                    continue

                dmat = distance.cdist(com1,com2)
                c_idx = (dmat<5.0*dd).nonzero()

                print len(c_idx[0])
                for idx in range(len(c_idx[0])):
                    # skip overlapping
                    #if(len(N.union1d(idx1[c_idx[0][idx]],idx2[c_idx[1][idx]]))!=l1+l2):
                    #    continue
                    ii = idx1[c_idx[0][idx]]+idx2[c_idx[1][idx]]
                    red_mat = mat[ii,:][:,ii]
                    ermsd = t.calc_dist_1d(ref_mat,red_mat)
                    if(ermsd < args.treshold):
                        print ermsd,
                        for el in ii:
                            print sequence[jj][el],
                        print ''

                    #grid0 = list(itertools.product(idx1[c_idx[0][idx]],idx2[c_idx[1][idx]]))
                    #grid1 = []
                    #for el in grid0:
                        
                    #grid1 = [x[0]*mat.shape[0]+x[1] for x in grid ]
                    #print grid0
                    #rr=  flat_mat[grid1]
                    

    if(args.type=='vector'):

        for ii,f in enumerate(files):
            atoms,sequence = pb.get_coord(f)

            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model."
                    print "# Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                assert(ref_mat.shape[0] == args.l1+args.l2)
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
                    mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                    for j in range(0,mat.shape[0]-l1+1):
                        ends =  N.union1d(mat[j,:].nonzero()[0],mat[:,j].nonzero()[0])
                        #ends =  N.intersection1d(mat[j,:].nonzero()[0],mat[:,j].nonzero()[0])
                        for k in ends:
                            # standard case
                            if(k<l2-1):
                                continue
                            start1 = j
                            end1 = j+l1-1
                            start2 = k-l2+1
                            end2 = k

                            if (all(mat[end1,start2]) == 0 and all(mat[start2,end1]) == 0):
                                continue
                            idx =  N.arange(start1,end1+1).tolist()
                            idx.extend(N.arange(start2,end2+1).tolist())
                            red_mat = mat[idx,:][:,idx]
                            #if(k<j):
                            #    red_mat = red_mat.T
                            ermsd = t.calc_dist_nd(ref_mat,red_mat)

                            if(ermsd < args.treshold):
                                seq = ' - '
                                for el in sequence[jj][idx[0]:idx[l1-1]+1]:
                                    seq += (el + ' ')
                                seq += " - "    
                                for el in sequence[jj][idx[l1]:idx[-1]+1]:
                                    seq += (el + ' ')
                                    
                                string = '%8.5f %s %i %s \n' % (ermsd,f,jj,seq)
                                fh.write(string)


                
####################### MOTIF ########################
