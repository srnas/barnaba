import mdtraj as md
import numpy as np
import sys
import nucleic
import definitions

def get_string(time,pairs,annotation,seq_id,hread=True):


    if(hread):
        string = "# t=%10s \n" % (time)
        ss = ["%10s %10s %2s \n" % (seq_id[pairs[j][0]],seq_id[pairs[j][1]],annotation[j]) for j in  xrange(len(pairs))]
        string += "".join(ss)
    else:
        #string = "%10s \n" % (time)
        string = ""
        #very expensive - use only for short sequences...
        merged = [str(el[0]) + "_" + str(el[1]) for el in pairs]
        for i1 in xrange(len(seq_id)):
            for i2 in xrange(i1+1,len(seq_id)):
                symbol = ".."
                lab = str(i1) + "_" + str(i2)
                if(lab in merged):
                    symbol = annotation[merged.index(lab)]
                string += "%3s " % symbol
        #string += ""
    return string

def wc_gaussian(vec):

     dev1  = vec - definitions.wc_mean
     maha1 = np.einsum('...k,...kl,...l->...', dev1, definitions.inv_sigma, dev1)
     return (2 * np.pi)**(-1.5)* definitions.det_sigma**(0.5)*np.exp(-0.5 * maha1)

def dihedral(p1,p2,p3,p4):

    # difference vectors b0 is reversed
    #b0 = vecs[:,1] - vecs[:,0]
    #b1 = vecs[:,1] - vecs[:,2]
    #b2 = vecs[:,2] - vecs[:,3]
    b0 = p2-p1
    b1 = p2-p3
    b2 = p3-p4
    # norm
    norm_sq = np.sum(b1**2)
    norm_sq_inv = 1.0/norm_sq

    #print (np.sum(b0*b1,axis=1)*b1).shape
    v0 = b0 - b1*((np.sum(b0*b1)*norm_sq_inv))
    v2 = b2 - b1*((np.sum(b0*b2)*norm_sq_inv))
    x = np.sum(v0*v2)
    m = np.cross(v0,b1)*np.sqrt(norm_sq_inv)
    y = np.sum(m*v2)
    return np.arctan2( y, x )

def dotbr(openings,closings,ll):

    # check pseudoknots
    dotbr = ['.']*ll
    levels = [-1]*ll
    for idx1 in xrange(len(openings)):
        start1 = openings[idx1]
        end1 = closings[idx1]
        # up one level
        levels[start1] += 1
        levels[end1] += 1
        
        for idx2 in xrange(len(closings)):
            end2 = closings[idx2]
            if(levels[end2] == levels[start1]):
                if(end2 > start1 and end2 < end1):
                    levels[start1] += 1
                    levels[end1] += 1

    for idx1 in xrange(len(openings)):
        start1 = openings[idx1]
        end1 = closings[idx1]
        
        dotbr[start1] = definitions.op[levels[start1]]
        dotbr[end1] = definitions.cl[levels[end1]]
    return "".join(dotbr)

def calculate(mat,angles,glyco,sequence):

    ll = len(sequence)
    # get upper triangular
    dist_sq = mat[:,:,0]**2 + mat[:,:,1]**2 + mat[:,:,2]**2
    # find pairs with cutoff < 1.6 (upper triangular only...)
    m_idx1 = np.where(np.triu(dist_sq,1)>0.0)
    m_idx2 = (m_idx1[1],m_idx1[0])
    rho1 = mat[m_idx1][:,0]**2 + mat[m_idx1][:,0]**2
    rho2 = mat[m_idx2][:,0]**2 + mat[m_idx2][:,0]**2
    z1 = mat[m_idx1][:,2]**2
    z2 = mat[m_idx2][:,2]**2
    
    openings = []
    closings = []
    pairs = []
    anno = []
    
    for i in range(len(m_idx1[0])):
        idx1=m_idx1[0][i]
        idx2=m_idx1[1][i]

        id1 = idx1,idx2
        id2 = idx2,idx1
        # stacking
        int_type = "XX"
        # criteria on |z| > 0.2

        if(z1[i]> 0.04 and z2[i] > 0.04):
            # criteria on angles (40 < x < 140)

            if(np.abs(angles[id1])>0.76):
                        
                # criteria on rho < 2.5 AA 
                if(rho1[i] < 0.0625 or rho2[i] < 0.0625): 
                    
                    if(mat[id1][2] > 0.02):
                        if(mat[id2][2] < -0.02):
                            int_type = ">>"
                        else:
                            int_type = "><"
                    else:
                        if(mat[id2][2] > 0.02):
                            int_type = "<<"
                        else:
                            int_type = "<>"
                            
        if(z1[i]< 0.04 and z2[i] < 0.04):
            int_type = ""
            edges = [np.arctan2(mat[id1][1],mat[id1][0]), np.arctan2(mat[id2][1],mat[id2][0])]
            for edge in edges:
                if(edge > definitions.theta1 and edge <= definitions.theta2):
                    int_type += "W"
                else:
                    if(edge > definitions.theta2 or edge <= definitions.theta3):
                        int_type += "H"
                    else:
                        int_type += "S"
            # watson-watson can be watson-crick
            if(int_type == "WW"):
                r1 = sequence[m_idx1[0][i]][0]
                r2 = sequence[m_idx1[1][i]][0]
                if(definitions.complementary[r1]==r2):
                    val = wc_gaussian(10.*mat[id1])*wc_gaussian(10.*mat[id2])
                    if(val>1.0e-08):
                        int_type = "WC"
                else:
                    if((r1=="G" and r2=="U") or (r1=="U" and r2=="G")):
                        int_type = "GU"

            # cis/trans
            ctrans = dihedral(glyco[1,idx1],glyco[0,idx1],glyco[0,idx2],glyco[1,idx2])
            if(np.abs(ctrans)<np.pi/2):
                int_type += "c"
            else:
                int_type += "t"

            # for generating dot-bracket
            if(int_type == "WCc"):
                assert (idx1 not in openings)
                assert (idx1 not in openings)
                assert (idx2 not in closings)
                assert (idx2 not in closings)
                openings.append(idx1)
                closings.append(idx2)
                
        pairs.append([idx1,idx2])
        anno.append(int_type)
    dotbracket = dotbr(openings,closings,ll)
    return pairs,anno, dotbracket

def annotate_traj(traj):
    top = traj.topology
    # initialize nucleic class
    nn = nucleic.Nucleic(top)
    rna_seq = ["%s_%s_%s" % (res.name,res.resSeq,res.chain.index) for res in nn.ok_residues]
    pairs = []
    annotations = []
    
    for i in xrange(traj.n_frames):

        coords_lcs = traj.xyz[i,nn.indeces_lcs]
        mat,angles = nn.get_mat_annotation(coords_lcs)
        if(len(mat)==0): continue
        
        glyco = traj.xyz[i,nn.indeces_glyco]
        pp, aa, dotbracket = calculate(mat,angles,glyco,nn.rna_seq_id)
        pairs.append(pp)
        annotations.append(aa)
        #stri = get_string(traj.time[i],pp,aa,rna_seq)
        #print stri
    return pairs, annotations

def annotate(filename,topology=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    pairs, annotations = annotate_traj(traj)
    return pairs, annotations
