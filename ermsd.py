import pdbreader as pb
import tools as t


def ermsd(args,files):
    

    fh = open(args.name,'w')

    print "# Calculating ERMSD..."
    pb.write_args(args,fh)

    if(args.type=='modulus'):

        # calculate interaction matrix for of the reference structure
        atoms,sequence = pb.get_coord(files[0])
        lcs_ref,origo_ref = t.coord2lcs(atoms[0])
        ref_mat = t.lcs2mat_1d(lcs_ref,origo_ref,args.cutoff)

        # all the rest - and calculate ERMSD on-the-fly
        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[0])

            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)
                assert(origo_ref.shape==origo.shape))               
                mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                ermsd = N.sqrt(sum(sum((ref_mat-mat)**2))/mat.shape[0])
                string = '%8.5f %s.%i \n' % (ermsd,f,jj)
                fh.write(string)
                

    if(args.type=='vector'):

        # calculate interaction matrix for of the reference structure
        atoms,sequence = pb.get_coord(files[0])
        lcs,origo = t.coord2lcs(atoms[0])
        ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)

        # all the rest - and calculate ERMSD on-the-fly
        for ii in xrange(1,len(files)):
            atoms,sequence = pb.get_coord(files[0])

            for jj,model in enumerate(atoms):
                lcs,origo = t.coord2lcs(model)

                assert(origo_ref.shape==origo.shape))               
                mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                diff = (mat-ref_mat)**2
                ermsd = N.sqrt(sum(sum(diff))/mat.shape[0])
                string = '%8.5f %s.%i \n' % (ermsd,f,jj)
                fh.write(string)
        

    fh.close()
    return 0


