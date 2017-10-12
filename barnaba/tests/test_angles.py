
import barnaba.barnaba as bb
import barnaba.definitions as dd

#fname = "data/1S72.pdb"
fname = "data/sample1.pdb"
traj = "data/samples.xtc"


angles_b,rr = bb.backbone_angles(fname)

fh = open("angles_0.dat",'w')
stri = "%20s " % "#"
for pp in dd.bb_angles:
    stri += " %10s " % pp
stri += "\n"
    
for e in range(angles_b.shape[1]):
    stri += "%20s " % (rr[e])
    for k in range(angles_b.shape[2]):
        stri += " %10.4f " % angles_b[0,e,k]
    stri += "\n"
fh.write(stri)
fh.close()


resi = ["RC5_1_0","RG_69_0"]
angles = ["alpha","beta","chi"]
angles_b,rr = bb.backbone_angles(fname,residue=resi,angle=angles)
stri = "%20s " % "#"
for pp in angles:
    stri += " %10s " % pp
stri += "\n"

for e in range(angles_b.shape[1]):
    stri += "%20s " % (rr[e])
    for k in range(angles_b.shape[2]):
        stri += " %10.4f " % angles_b[0,e,k]
    stri += "\n"
    
fh = open("angles_1.dat",'w')
fh.write(stri)
fh.close()


resi = ["RG_69_0","RU_37_0"]
angles = ["gamma"]
angles_b,rr = bb.backbone_angles(traj,topology=fname,residue=resi,angle=angles)
stri = ""
for p in range(angles_b.shape[0]):
    for k in range(angles_b.shape[2]):
        stri += " %10.4f %10.4f " % (angles_b[p,0,k],angles_b[p,1,k])
    stri += "\n"
    
fh = open("angles_2.dat",'w')
fh.write(stri)
fh.close()


