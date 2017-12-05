import barnaba as bb

fname = "data/UUCG.xtc"
tname = "data/UUCG.pdb"

stackings, pairings, res = bb.annotate(fname,tname)
dotbr = bb.dot_bracket(pairings,res)

fh = open("dotbracket.test.dat",'w')
stri = ""
for k in range(len(pairings)):
    stri += " %06d " % k
    stri += "%s \n" % dotbr[k]
fh.write(stri)

fh.close()








