import sys


print "# This is a helper script"
print "# USAGE: python motif_helper.py filename.MOTIF.rna"

# read barnaba output first
barnaba = []
fh = open(sys.argv[1])
for line in fh:
    if("#" not in line):
         lista = []
         for ll in line.split("-")[1:]:
             lista.extend(ll.split())
         d = [line.split()[1],int(line.split()[2]),lista]
         barnaba.append(d)
fh.close()

root = sys.argv[1].split(".rna")[0]

for ii,el in enumerate(barnaba):

    pdb_out = root + "_" + str(ii) + ".pdb"
    pdb = open(el[0],'r')
    fh_out = open(pdb_out,'w')
    string = "REMARK " + el[0] + '\n'
    string += "REMARK " + str(el[2]) + "\n"
    print ""
    fh_out.write(string)

    for line in pdb:

        if(line.split()[0] != "ATOM"):
            continue
        
        RESN = line[22:26].strip()
        CHAIN = line[21:22].strip()
        REST = line[17:20].strip()
        RES_ID = RESN + "_" + REST + "_" + CHAIN
        if(RES_ID in el[2]):
            fh_out.write(line)

    pdb.close()
    fh_out.close()
