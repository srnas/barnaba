
import barnaba.barnaba as bb

#fname = "data/sample1.pdb"
fname ="/Users/sandrobottaro/Projects/OPC/notebook/Tetranucleotides/AAAA/AAAAh.pdb"
traj="/Users/sandrobottaro/Projects/OPC/notebook/Tetranucleotides/reweight/AAAA_TIP3P_300/traj_AAAA_TIP3P_300.dcd"
angles,rr = bb.jcouplings(traj,topology=fname)
#angles,rr = bb.jcouplings(fname,couplings=["H3P","H2H3","H1H2"])
#angles,rr = bb.jcouplings(fname)
print angles.shape
fh = open("jcouplings.dat",'w')
for e in range(angles.shape[0]):
    stri = ""
    for k in range(angles.shape[1]):
        for l in range(angles.shape[2]):
            stri += " %10.4f " % angles[e,k,l]
        stri += "\n"
    stri += "\n"
    fh.write(stri)
fh.close()
#exit()

angles,rr = bb.jcouplings(traj,topology=fname,raw=True)
#angles,rr = bb.jcouplings(fname,raw=True,couplings=["H3P","H2H3","H1H2"])
#angles,rr = bb.jcouplings(fname,raw=True)
print angles.shape
fh = open("jcouplings_raw.dat",'w')
for e in range(angles.shape[0]):
    stri = ""
    for k in range(angles.shape[1]):
        for l in range(angles.shape[2]):
            stri += " %10.4f " % angles[e,k,l]
        stri += "\n"
    stri += "\n"
    fh.write(stri)
fh.close()
