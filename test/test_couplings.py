import barnaba as bb

top = "data/sample1.pdb"
traj = "data/samples.xtc"


angles,rr = bb.jcouplings(traj,topology=top)
fh = open("jcouplings.test.dat",'w')
for e in range(angles.shape[0]):
    stri = ""
    for k in range(angles.shape[1]):
        for l in range(angles.shape[2]):
            stri += " %10.4f " % angles[e,k,l]
        stri += "\n"
    stri += "\n"
    fh.write(stri)
fh.close()

angles,rr = bb.jcouplings(traj,topology=top,raw=True)
fh = open("jcouplings_raw.test.dat",'w')
for e in range(angles.shape[0]):
    stri = ""
    for k in range(angles.shape[1]):
        for l in range(angles.shape[2]):
            stri += " %10.4f " % angles[e,k,l]
        stri += "\n"
    stri += "\n"
    fh.write(stri)
fh.close()
