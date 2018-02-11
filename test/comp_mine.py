from __future__ import absolute_import, division, print_function

def comp(filename):

    def readfile(ff):
        data = []
        fh = open(ff)
        for line in fh:
            ss = line.split()
            if("REMARK" in line): continue
            if(len(ss)!=0):
                vv = []
                for el in ss:
                    if(el=="nan"):
                        vv.append(el)
                    else:
                        try:
                            vv.append(float(el))
                        except:
                            vv.append(el)
                data.append(vv)
        fh.close()
        return data
    
    ref = readfile(filename)
    cur = readfile(filename.replace("reference","tmp"))
    assert len(ref)==len(cur), "reference and test files %s have different shapes " % filename.split("/")[-1]
    for el1,el2 in zip(ref,cur):
        assert len(el1)==len(el2), "reference and test files %s have different shapes " % filename.split("/")[-1]
        for it1,it2 in zip(el1,el2):
            if(type(it1)==float):
                assert (it1-it2)**2<1.E-05, "reference and test files %s give different results " % filename.split("/")[-1]
            else:
                assert it1==it2, "reference and test files %s give different results " % filename.split("/")[-1]
