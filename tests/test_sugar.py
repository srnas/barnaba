import barnaba as bb
import barnaba.definitions as dd

fname = "data/sample1.pdb"
traj = "data/samples.xtc"

angles,rr = bb.sugar_angles(fname)
fh = open("sugar_0.test.dat",'w')
stri = "%20s " % "#"
for pp in dd.sugar_angles:
    stri += " %10s " % pp
stri += "\n"
    
for e in range(angles.shape[1]):
    stri += "%20s " % (rr[e])
    for k in range(angles.shape[2]):
        stri += " %10.4f " % angles[0,e,k]
    stri += "\n"
fh.write(stri)
fh.close()


resi = ["RC5_1_0","RG_69_0","RU_37_0"]
angles = ["nu1","nu4","nu3"]
angles_b,rr = bb.sugar_angles(fname,residues=resi,angles=angles)
stri = "%20s " % "#"
for pp in angles:
    stri += " %10s " % pp
stri += "\n"

for e in range(angles_b.shape[1]):
    stri += "%20s " % (rr[e])
    for k in range(angles_b.shape[2]):
        stri += " %10.4f " % angles_b[0,e,k]
    stri += "\n"
    
fh = open("sugar_1.test.dat",'w')
fh.write(stri)
fh.close()

resi = ["RG_69_0","RU_37_0"]
angles = ["nu4","nu2"]
angles_b,rr = bb.sugar_angles(traj,topology=fname,residues=resi,angles=angles)
stri = ""
for p in range(angles_b.shape[0]):
    for k in range(angles_b.shape[2]):
        stri += " %10.4f %10.4f " % (angles_b[p,0,k],angles_b[p,1,k])
    stri += "\n"
    
fh = open("sugar_2.test.dat",'w')
fh.write(stri)
fh.close()



angles,rr = bb.pucker_angles(fname)

fh = open("sugar_3.test.dat",'w')
stri = "%20s " % "#"
na = ["Phase","tm"]
for pp in na:
    stri += " %10s " % pp
stri += "\n"
    
for e in range(angles.shape[1]):
    stri += "%20s " % (rr[e])
    for k in range(angles.shape[2]):
        stri += " %10.4f " % angles[0,e,k]
    stri += "\n"
fh.write(stri)
fh.close()


resi = ["RC5_1_0","RG_69_0","RU_37_0"]
angles_b,rr = bb.pucker_angles(fname,residues=resi)
stri = "%20s " % "#"
for pp in na:
    stri += " %10s " % pp
stri += "\n"

for e in range(angles_b.shape[1]):
    stri += "%20s " % (rr[e])
    for k in range(angles_b.shape[2]):
        stri += " %10.4f " % angles_b[0,e,k]
    stri += "\n"
    
fh = open("sugar_4.test.dat",'w')
fh.write(stri)
fh.close()
