from __future__ import absolute_import, division, print_function
import sys

# read new
new_anno = []
new_anno1 = []

fh = open(sys.argv[1])
for line in fh:
    if("#" not in line):
        new_anno.append(line.split()[0] + "/" + line.split()[1])
        new_anno1.append(line.split()[2])
fh.close()

# read old
old_anno = []
old_anno1 = []
fh = open(sys.argv[2])
for line in fh:
    if("#" not in line):
        old_anno.append(line.split()[0] + "/" + line.split()[1])
        old_anno1.append(line.split()[2])
fh.close()

oks =0
for j in range(len(new_anno)):
    if(new_anno[j] in old_anno):
        k = old_anno.index(new_anno[j])
        if(old_anno1[k]==new_anno1[j]):
            oks += 1
        else:
            print "MMH %s has different annotation:new %s old %s " % (new_anno[j], new_anno1[j],old_anno1[k])
    else:
        print "MMH %s not present in old %s " % (new_anno[j], new_anno1[j])

print "oks: %d (%d)" % (oks,len(new_anno))
print "#######"


oks =0
for j in range(len(old_anno)):
    if(old_anno[j] in new_anno):
        k = new_anno.index(old_anno[j])
        if(new_anno1[k]==old_anno1[j]):
            oks += 1
        else:
            print "MMH %s has different annotation:old %s new %s " % (old_anno[j], old_anno1[j],new_anno1[k])
    else:
        print "MMH %s not present in new" % oldanno[j]

print "oks: %d (%d)" % (oks,len(old_anno))
print "#######"
