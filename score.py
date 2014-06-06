import pdbreader as pb
import kde as kde
import tools as t

####################### SCORING ########################
        
def score(args,files):

    fh = open(args.name,'w')

    pb.write_args(args,fh)

    ref_atoms,ref_sequence = pb.get_coord(args.ff)
    ref_lcs,ref_origo = t.coord2lcs(ref_atoms[0])
    ref_mat,ids = t.lcs2mat_score(ref_lcs,ref_origo,args.cutoff)
    
    kernel = kde.gaussian_kde(ref_mat)
    kernel.set_bandwidth(0.25)

    print "# KDE computed. Bandwidth=",kernel.factor
    print "# Calculating ESCORE..."

    for f in files:
        atoms,sequence = pb.get_coord(f)
        for ii,model in enumerate(atoms):
            lcs,origo = t.coord2lcs(model)

            # Cutoff is slightly augmented 
            mat,ids = t.lcs2mat_score(lcs,origo,args.cutoff+0.2)
            val = kernel(mat)
            #for kk in range(len(sequence[ii])):
            #    ss = 0.0
            #    for el in range(len(val)):
            #        if(kk in ids[el]):
            #            ss += val[el]
            #    print sequence[ii][kk],ss

            string = '%8.5f %s.%i \n' % (sum(val),f,ii)
            fh.write(string)

####################### SCORING ########################
