

mmap = {"0":"0","1":"9"}

def compare(d1,a1,d2,a2):
    
    common = 0
    common_anno = 0
    d1_only = 0
    d2_only = 0
    
    for j in range(len(d1)):
        if(d1[j] in d2):
            k = d2.index(d1[j])
            common_anno += 1
            if(a2[k]==a1[j]):
                common += 1
        else:
            d1_only += 1
    common1 = 0
    common_anno1 = 0
                
    for j in range(len(d2)):
        if(d2[j] in d1):
            k = d1.index(d2[j])
            common_anno1 += 1
            if(a1[k]==a2[j]):
                common1 += 1
        else:
            d2_only += 1
    assert(common1==common)
    assert(common_anno1==common_anno)
    return common, common_anno, d1_only, d2_only
    
# barnaba annotation ##
fh = open("pairings.dat")
b_anno = []
b_anno1 = []
for line in fh:
    if("#" not in line):
        ll1 = line.split()[0]
        ll2 = line.split()[1]
        r1 = "%s.%s%s" % (mmap[ll1.split("_")[2]],ll1.split("_")[0],ll1.split("_")[1])
        r2 = "%s.%s%s" % (mmap[ll2.split("_")[2]],ll2.split("_")[0],ll2.split("_")[1])
        rr = "/".join(sorted([r1,r2]))
        anno = line.split()[2]
        if(anno=="XXX"): continue
        if(anno=="WCc" or anno=="GUc"):
            anno = "WWc"
        b_anno1.append(anno)
        b_anno.append(rr)
fh.close()
print "# read %d barnaba pairing" % len(b_anno)

# read dssr annotation 
d_anno = []
d_anno1 = []
fh = open("1S72_dssr.anno")
for line in fh:
    r1 = line.split()[1]
    r2 = line.split()[2]
    ll = "/".join(sorted([r1,r2]))
    #if(ll in d_anno): continue
    anno = line.split()[6]
    anno1 = anno[1:] + anno[0]
    d_anno1.append(anno1)
    d_anno.append(ll)    
fh.close()
print "# read %d dssr pairing" % len(d_anno)

# read fred annotation
f_anno = []
f_anno1 = []
fh = open("1S72_fr3d_basepairs.csv")
for line in fh:
    r1 = line.split(",")[0][1:-1].split("|")
    r2 = line.split(",")[2][1:-2].split("|")
    rr1 = "%s.%s%s" % (r1[2],r1[3],r1[4])
    rr2 = "%s.%s%s" % (r2[2],r2[3],r2[4])
    ll = "/".join(sorted([rr1,rr2]))
    if(ll in f_anno): continue
    anno = line.split(",")[1][1:-1]
    if(anno[0] == "n"):
        anno = anno[1:]
    anno1 = anno[1:] + anno[0] 
    f_anno1.append(anno1)
    f_anno.append(ll)
print "# read %d fr3d pairing" % len(f_anno)


print "Fred-dssr"
print compare(f_anno,f_anno1,d_anno,d_anno1)

print "Fred-barnaba"
print compare(f_anno,f_anno1,b_anno,b_anno1)

print "barnaba-dssr"
print compare(b_anno,b_anno1,d_anno,d_anno1)
