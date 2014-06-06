from scipy.spatial import distance
import numpy as N
import pdbreader as pb
import tools as t


####################### MOTIF #########################
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

def ss_motif(args,files):

    fh = open(args.name,'w')

    print "# Finding 3D motifs..."
    pb.write_args(args,fh)

    atoms,sequence = pb.get_coord(files[0])
    lcs,origo = t.coord2lcs(atoms[0])
    ll = len(lcs)

    if(args.type=='modulus'):
        
        ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
        ref_mat = ref_mat.reshape(-1)
        
        # loop over structures
        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])
            
            # loop over models
            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                mat = t.lcs2mat_1d(lcs,origo,args.cutoff)

                idx = []
                # loop in submatrices
                for j in range(0,mat.shape[0]-ll+1):
                    tmp_idxs = do_idx(j,ll,args.bulges)
                    for tmp_idx in tmp_idxs:
                        if(tmp_idx[-1]>=mat.shape[0] or tmp_idx in idx):
                            continue
                        
                        red_mat = mat[tmp_idx,:][:,tmp_idx]
                        red_mat = red_mat.reshape(-1)
                        ermsd = N.sqrt( sum((ref_mat-red_mat)**2)/ll)
                        idx.append(tmp_idx)
                        if(ermsd < args.treshold):
                            seq = ' - '
                            for el in tmp_idx:
                                seq += (sequence[jj][el] + ' ')
                            string = '%8.5f %s %i %s \n' % (ermsd,files[ii],jj,seq)
                            fh.write(string)


    if(args.type=='vector'):

        ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
        ref_mat = ref_mat.reshape(-1,4)

        # loop over structures
        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])

            # loop over models
            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                idx = []
                # loop in submatrices
                for j in range(0,mat.shape[0]-ll+1):
                    tmp_idxs = do_idx(j,ll,args.bulges)
                    for tmp_idx in tmp_idxs:
                        if(tmp_idx[-1]>=mat.shape[0] or tmp_idx in idx):
                            continue
                        
                        red_mat = mat[tmp_idx,:][:,tmp_idx]
                        red_mat = red_mat.reshape(-1,4)
                        diff = (red_mat-ref_mat)**2
                        ermsd = N.sqrt(sum(sum(diff))/ll)
                        idx.append(tmp_idx)
                        if(ermsd < args.treshold):
                            seq = ' - '
                            for el in tmp_idx:
                                seq += (sequence[jj][el] + ' ')
                            string = '%8.5f %s %i %s \n' % (ermsd,files[ii],jj,seq)
                            fh.write(string)

