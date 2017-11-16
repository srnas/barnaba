import numpy as np
import barnaba.barnaba as bb
import sys

aa, rr = bb.jcouplings(sys.argv[1],raw=True,couplings=['H1H2', 'H2H3', 'H3H4'])
ff=180./np.pi
for i in range(len(rr)):
    #print rr[i],np.average(aa[0,:,i]),np.average(aa[1,:,i]),np.average(aa[2,:,i])
    print rr[i],ff*aa[0,0,i],ff*aa[1,0,i],ff*aa[2,0,i]
