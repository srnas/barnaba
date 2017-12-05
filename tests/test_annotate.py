
import barnaba.barnaba as bb


fname = "data/samples.xtc"
tname = "data/sample1.pdb"
stackings, pairings, res = bb.annotate(fname,tname)

fh = open("stackings_traj.dat",'w')
stri = "# STACKING \n"
for k in range(len(stackings)):
    stri += "# frame %d \n" % k
    for e in range(len(stackings[k][0])):
        stri += "%15s " % (res[stackings[k][0][e][0]])
        stri += "%15s " % (res[stackings[k][0][e][1]])
        stri += " %4s \n" % (stackings[k][1][e])
    stri += "\n"
fh.write(stri)
fh.close()

fh = open("pairings_traj.dat",'w')
stri = "# PAIRING \n"
for k in range(len(pairings)):
    stri += "# frame %d \n" % k
    for e in range(len(pairings[k][0])):
        stri += "%15s " % (res[pairings[k][0][e][0]])
        stri += "%15s " % (res[pairings[k][0][e][1]])
        stri += " %4s \n" % (pairings[k][1][e])
fh.write(stri)
fh.close()



fname = "data/1S72.pdb"
stackings, pairings, res = bb.annotate(fname)
fh = open("stackings.dat",'w')
stri = "# STACKING \n"
for e in range(len(stackings[0][0])):
    stri += "%15s " % (res[stackings[0][0][e][0]])
    stri += "%15s " % (res[stackings[0][0][e][1]])
    stri += " %4s \n" % (stackings[0][1][e])
fh.write(stri)
fh.close()

fh = open("pairings.dat",'w')
stri = "# PAIRING \n"
for e in range(len(pairings[0][0])):
    stri += "%15s " % (res[pairings[0][0][e][0]])
    stri += "%15s " % (res[pairings[0][0][e][1]])
    stri += " %4s \n" % (pairings[0][1][e])
fh.write(stri)
fh.close()




