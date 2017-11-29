
import barnaba.barnaba as bb


fname = "data/samples.xtc"
tname = "data/sample1.pdb"
fname = "data/UUCG.xtc"
tname = "data/UUCG.pdb"

stackings, pairings, res = bb.annotate(fname,tname)
dotbr = bb.dot_bracket(pairings,res)

fh = open("dotbracket.dat",'w')
stri = ""
for k in range(len(pairings)):
    stri += " %06d " % k
    #for e in range(len(pairings[k][0])):
    #    if(pairings[k][1][e]=="WCc" or pairings[k][1][e]=="GUc"):
    #        stri += "%15s " % (res[pairings[k][0][e][0]])
    #        stri += "%15s " % (res[pairings[k][0][e][1]])
    #        stri += " %4s \n" % (pairings[k][1][e])
    stri += "%s \n" % dotbr[k]
fh.write(stri)

fh.close()








