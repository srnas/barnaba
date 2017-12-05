import barnaba as bb

fname = "data/sample1.pdb"
rvecs,resi = bb.dump_rvec(fname,cutoff=1.7)

fh = open("dump_r.test.dat",'w')
stri = ""
for i1 in range(len(resi)):
    for i2 in range(len(resi)):
        if(sum((rvecs[0,i1,i2])**2) > 1.0e-5):
            stri += "%10s %10s %14e %14e %14e \n" % (resi[i1],resi[i2],rvecs[0,i1,i2,0],rvecs[0,i1,i2,1],rvecs[0,i1,i2,2])
fh.write(stri)
fh.close()


gvecs,resi = bb.dump_gvec(fname)
fh = open("dump_g.test.dat",'w')
stri = ""
for i1 in range(len(resi)):
    for i2 in range(len(resi)):
        if(sum((gvecs[0,i1,i2])**2) > 1.0e-5):
            stri += "%10s %10s %14e %14e %14e %14e  \n" % (resi[i1],resi[i2],gvecs[0,i1,i2,0],gvecs[0,i1,i2,1],gvecs[0,i1,i2,2],gvecs[0,i1,i2,3])
fh.write(stri)
fh.close()
