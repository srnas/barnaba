from __future__ import absolute_import, division, print_function
import sys
mmap = {"0":"0","9":"1"}
# read new
new_anno = []
new_anno1 = []

fh = open(sys.argv[1])
for line in fh:
    if("#" not in line):
        rr = "/".join(sorted([line.split()[0],line.split()[1]]))
        new_anno.append(rr)
        anno = line.split()[2]
        if(anno=="XXX"):continue
        if(anno=="WCc" or anno=="GUc"):
            anno = "WWc"
        new_anno1.append(anno)
fh.close()

# read csv
old_anno = []
old_anno1 = []
fh = open(sys.argv[2])
for line in fh:
    r1 = line.split(",")[0][1:-1]
    r2 = line.split(",")[2][1:-2]
    rr1 = "%s_%s_%s" % (r1.split("|")[3],r1.split("|")[4],mmap[r1.split("|")[2]])
    rr2 = "%s_%s_%s" % (r2.split("|")[3],r2.split("|")[4],mmap[r2.split("|")[2]])
    ll = "/".join(sorted([rr1,rr2]))
    if(ll in old_anno): continue
    anno = line.split(",")[1][1:-1]
    if(anno[0] == "n"):
        anno = anno[1:]
    anno1 = anno[1:] + anno[0] 
    old_anno1.append(anno1)
    old_anno.append(ll)
    
fh.close()

oks =0
for j in range(len(new_anno)):
    if(new_anno[j] in old_anno):
        k = old_anno.index(new_anno[j])
        if(old_anno1[k]==new_anno1[j]):
            oks += 1
        #else:
        #    print "MMH %s has different annotation:new %s old %s " % (new_anno[j], new_anno1[j],old_anno1[k])
    else:
        print "MMH %s not present in fred %s " % (new_anno[j], new_anno1[j])

print "oks: %d (%d)" % (oks,len(new_anno))
print "#######"


oks =0
for j in range(len(old_anno)):
    if(old_anno[j] in new_anno):
        k = new_anno.index(old_anno[j])
        if(new_anno1[k]==old_anno1[j]):
            oks += 1
        #else:
        #    print "MMH %s has different annotation:old %s new %s " % (old_anno[j], old_anno1[j],new_anno1[k])
    else:
        print "MMH %s not present in mine %s " % (old_anno[j], old_anno1[j])

print "oks: %d (%d)" % (oks,len(old_anno))
print "#######"
