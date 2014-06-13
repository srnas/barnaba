import pdbreader as pb
import tools as t
import numpy as N

def ermsd(args,files):
    

    fh = open(args.name,'w')

    print "# Calculating ERMSD..."
    pb.write_args(args,fh)

    # calculate interaction matrix of the reference structure
    atoms,sequence = pb.get_coord(files[0])
    lcs_ref,origo_ref = t.coord2lcs(atoms[0])

    if(args.type=='modulus'):

        # calculate interaction matrix of the reference structure
        ref_mat = t.lcs2mat_1d(lcs_ref,origo_ref,args.cutoff)
        ref_mat_f = ref_mat.reshape(-1)

        # all the rest - and calculate ERMSD on-the-fly
        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])

            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                assert(origo_ref.shape==origo.shape)               
                mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                mat_f = mat.reshape(-1)
                ermsd = N.sqrt( sum((ref_mat_f-mat_f)**2)/len(lcs))
                if(args.ermsf==False):
                    string = '%8.5f %s.%i \n' % (ermsd,files[ii],jj)
                else:
                    string = '%8.5f - ' % (ermsd)
                    for k in xrange(len(lcs)):
                        ermsf = 0.5*(N.sqrt( sum((ref_mat[k,:]-mat[k,:])**2)/len(lcs)) + N.sqrt( sum((ref_mat[:,k]-mat[:,k])**2)/len(lcs)))
                        string += " %8.5f " % (ermsf)
                    string += '- %s.%i \n' % (files[ii],jj)
                fh.write(string)
                

    if(args.type=='vector'):

        ref_mat = t.lcs2mat_4d(lcs_ref,origo_ref,args.cutoff)
        ref_mat_f = ref_mat.reshape(-1,4)
        
        # all the rest - and calculate ERMSD on-the-fly
        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[ii])
            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)

                assert(origo_ref.shape==origo.shape)               
                mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                mat_f = mat.reshape(-1,4)
                diff = (mat_f-ref_mat_f)**2
                ermsd = N.sqrt(sum(sum(diff))/len(lcs))
                if(args.ermsf==False):
                    string = '%8.5f %s.%i \n' % (ermsd,files[ii],jj)
                else:
                    string = '%8.5f - ' % (ermsd)
                    for k in xrange(len(lcs)):
                        diff1 = (mat[k,:]-ref_mat[k,:])**2
                        diff2 = (mat[:,k]-ref_mat[:,k])**2
                        ermsf = 0.5*(N.sqrt( sum(sum(diff1))/len(lcs)) + N.sqrt( sum(sum(diff2))/len(lcs)))
                        string += " %8.5f " % (ermsf)
                    string += '- %s.%i \n' % (files[ii],jj)
                fh.write(string)
        

    fh.close()
    return 0


