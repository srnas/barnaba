
import barnaba.angles as bb

fname = "data/1S72.pdb"

rr, angles_s = bb.sugar_angles(fname)
rr, angles_b = bb.backbone_angles(fname)

fh = open("angles.dat",'w')
stri = ""
for e in range(len(rr)):
    stri += "%20s " % (rr[e])
    for k in range(7):
        stri += " %10.4f " % angles_b[0,e,k]
    for k in range(7):
        stri += " %10.4f " % angles_s[0,e,k]
    stri += "\n"
fh.write(stri)
fh.close()
