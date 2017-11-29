import barnaba.barnaba as bb
import sys
import mdtraj as md

tname = "data/UUCG.pdb"
fname = "data/UUCG.xtc"
t = md.load(fname, top=tname)
t[2330].save_pdb("crap.pdb")
exit()
rvecs,resi = bb.dump_rvec(fname,tname,cutoff=1.7)


aa = [2331,2330,2332]
for j in aa:
    print j
    for i1 in range(len(resi)):
        for i2 in range(len(resi)):
            if((resi[i1]=="C_2_0" and resi[i2] == "G_6_0") and sum((rvecs[j,i1,i2])**2)>0.0001):
                print "%10s %10s %14e %14e %14e \n" % (resi[i1],resi[i2],rvecs[j,i1,i2,0],rvecs[j,i1,i2,1],rvecs[j,i1,i2,2])
            if((resi[i1]=="G_6_0" and resi[i2] == "C_2_0") and sum((rvecs[j,i1,i2])**2)>0.0001):
                print "%10s %10s %14e %14e %14e \n" % (resi[i1],resi[i2],rvecs[j,i1,i2,0],rvecs[j,i1,i2,1],rvecs[j,i1,i2,2])
            if((resi[i1]=="C_2_0" and resi[i2] == "G_7_0") and sum((rvecs[j,i1,i2])**2)>0.0001):
                print "%10s %10s %14e %14e %14e \n" % (resi[i1],resi[i2],rvecs[j,i1,i2,0],rvecs[j,i1,i2,1],rvecs[j,i1,i2,2])
            if((resi[i1]=="G_7_0" and resi[i2] == "C_2_0") and sum((rvecs[j,i1,i2])**2)>0.0001):
                print "%10s %10s %14e %14e %14e \n" % (resi[i1],resi[i2],rvecs[j,i1,i2,0],rvecs[j,i1,i2,1],rvecs[j,i1,i2,2])
            

