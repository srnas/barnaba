import pdbreader as pb
import tools as t
import numpy as N

##################### ANNOTATE #######################

def annotate(args,files):

    print "# Annotating RNA structures..."

    fh = open(args.name,'w')
    pb.write_args(args,fh)

    for f in files:
        atoms,sequence = pb.get_coord(f)
        for model,seq in zip(atoms,sequence):
            lcs,origo = t.coord2lcs(model)
            mat = t.lcs2mat_3d(lcs,origo,args.cutoff)
            int_mat = t.analyze_mat(mat,seq)
            #idx = N.triu_indices(len(seq),1)
            #data.append(int_mat[idx])
            
    #data = N.array(data)
                # remove columns with no interactions and print to file
    #idx1 = []
    #header = "# "
    #for k in range(data.shape[1]):
    #    if(any(data[:,k]) == 0):
    #        continue
    #    else:
    #        header += ref_sequence[0][idx[0][k]] + "/" + ref_sequence[0][idx[1][k]] + ":" + str(len(idx1)+2) + " "
    #        idx1.append(k)
    #header += "\n"
    #fh.write(header)
    oo = [1.0,2.0,3.0,4.0,5.0]
    # write data to file
    for i in range(mat.shape[0]/2):
        m1 = int_mat[i,:].tolist()
        if(5.0 in m1):
            idx = m1.index(5.0)
            if(int_mat[i+1,idx-1]==5.0):
                
                
                # WC
#                print "%15s %15s " % (seq[i],seq[idx]),
#                for ee in mat[i,idx]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[i,idx])]#
#
#                print "%15s %15s " % (seq[idx],seq[i]),
#                for ee in mat[idx,i]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[i,idx])]
#
#
#                print "%15s %15s " % (seq[i+1],seq[idx-1]),
#                for ee in mat[i+1,idx-1]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[i+1,idx-1])]
#
#
#                print "%15s %15s " % (seq[idx-1],seq[i+1]),
#                for ee in mat[idx-1,i+1]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[i+1,idx-1])]
#                ###
                
                
                # stacking 5->3
#                print "%15s %15s " % (seq[i],seq[i+1]),
#                for ee in mat[i,i+1]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[i,i+1])]
#
#                print "%15s %15s " % (seq[idx-1],seq[idx]),
#                for ee in mat[idx-1,idx]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[idx-1,idx])]

                
#                # stacking 3->5
#                print "%15s %15s " % (seq[i+1],seq[i]),
#                for ee in mat[i+1,i]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[i,i+1])]
#
#
#                print "%15s %15s " % (seq[idx],seq[idx-1]),
#                for ee in mat[idx,idx-1]:
#                    print "%8.4f " % ee,
#                print t.interactions[int(int_mat[idx-1,idx])]
#

                # mancano i due a croce
                #print "%15s %15s " % (seq[i],seq[idx-1]),
                #for ee in mat[i,idx-1]:
                #    print "%8.4f " % ee,
                #print t.interactions[int(int_mat[i,idx-1])]
                #continue
                
                #print "%15s %15s " % (seq[idx-1],seq[i]),
                #for ee in mat[idx-1,i]:
                #    print "%8.4f " % ee,
                #print t.interactions[int(int_mat[i,idx-1])]
                
               

                #print "%15s %15s " % (seq[i+1],seq[idx]),
                #for ee in mat[i+1,idx]:
                #    print "%8.4f " % ee,
                #print t.interactions[int(int_mat[i+1,idx])]

                
                print "%15s %15s " % (seq[idx],seq[i+1]),
                for ee in mat[idx,i+1]:
                    print "%8.4f " % ee,
                print t.interactions[int(int_mat[i+1,idx])]
                continue



    #for k in range(data.shape[0]):
    #    
    #    string = "%6i" % k
    ##    for l in idx1:
    #        string += "%4s " % (t.interactions[int(data[k,l])])
    #    print string
    #            if(all(mat[i,k]) == 0):
     #                   continue
     #               print seq[i].split("_")[1]+seq[k].split("_")[1],

            


    fh.close()
    return 0
            
##################### ANNOTATE #######################
