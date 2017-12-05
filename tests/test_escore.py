import barnaba.escore as escore

fname = "data/1S72.pdb"
# set "force-field"
ee = escore.Escore([fname])


fname = "data/samples.xtc"
tname = "data/sample1.pdb"
ss =  ee.score(fname,tname)
fh = open("escore.test.dat",'w')
stri = ""
for e in range(len(ss)):
    stri += "%10.4e \n " % (ss[e])
fh.write(stri)
fh.close()
