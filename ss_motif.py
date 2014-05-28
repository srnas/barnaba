from scipy.spatial import distance
import numpy as N
import pdbreader as pb
import tools as t


####################### MOTIF #########################

def ss_motif(args,files):

    fh = open(args.name,'w')

    print "# Finding 3D motifs..."
    pb.write_args(args,fh)

    
    if(args.type=='modulus'):

        for ii,f in enumerate(files):

            atoms,sequence = pb.get_coord(f)
            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model."
                    print "# Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                ll = ref_mat.shape[0]
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
      
                    mat = t.lcs2mat_1d(lcs,origo,args.cutoff)
                    for j in range(0,mat.shape[0]-ll+1):
                        red_mat = mat[j:j+ll,j:j+ll]
                        ermsd = t.calc_dist_1d(ref_mat,red_mat)
                        if(ermsd < args.treshold):
                            seq = ' - '
                            for el in sequence[jj][j:j+ll]:
                                seq += (el + ' ')
                            string = '%8.5f %s %i %s \n' % (ermsd,f,jj,seq)
                            fh.write(string)


    if(args.type=='vector'):

        for ii,f in enumerate(files):

            atoms,sequence = pb.get_coord(f)
            if(ii==0):
                if(len(atoms)!=1):
                    print "# Warning: your reference PDB file contains more than one model."
                    print "# Only the first one will be considered"
                lcs,origo = t.coord2lcs(atoms[0])
                ref_mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                ll = ref_mat.shape[0]
            else:
                for jj,model in enumerate(atoms):
                    lcs,origo = t.coord2lcs(model)
      
                    mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                    for j in range(0,mat.shape[0]-ll+1):
                        red_mat = mat[j:j+ll,j:j+ll]
                        ermsd = t.calc_dist_nd(ref_mat,red_mat)
                        if(ermsd < args.treshold):
                            string = '%8.5f %s %i %s \n' % (ermsd,f,jj,sequence[jj][j:j+ll])
                            fh.write(string)

