#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2014 Sandro Bottaro (sbottaro@sissa.it)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pdbreader as pb
import tools as t
import  numpy as N

def stringify(fname,jj,data,seq,hread=True):

    if(hread==False):
        data = data.reshape(-1)
        s = ''
        for el in data:
            s += '%10.6f ' % el
        s += '%s.%i \n' % (fname,jj)
        return s
    else:
        print len(seq),data.shape
        s = '# %s %i \n' % (fname,jj)
        for j in range(data.shape[0]):
            for k in range(data.shape[0]):
                if(j!=k and all(data[j,k]) != 0.0):
                    s += '%10s %10s ' % (seq[j],seq[k])
                    for l in range(len(data[j,k])):
                        s += '%10.6f ' % data[j,k,l]
                    s+= "\n"
        return s
    
def dump(args,files):
    

    print "# Calculating DUMP ... (ahah)"

    if(args.dumpG==True):
        fh_dumpG = open(args.name+ ".gvec",'w')

    if(args.dumpR==True):
        fh_dumpR = open(args.name+ ".rvec",'w')

    if(args.dumpL==True):
        fh_dumpL = open(args.name+ ".lcs",'w')
        
    for ii in xrange(0,len(files)):
        atoms,sequence = pb.get_coord(files[ii])
        
        for jj,model in enumerate(atoms):
            lcs,origo = t.coord2lcs(model)
            if(args.dumpL==True):

                s = '# %s %i \n' % (files[ii],jj)
                for j in range(lcs.shape[0]):
                    s += '%10s ' % sequence[jj][j]
                    for k in range(0,3):
                        for l in range(3):
                            s += '%10.6f ' % lcs[j,k,l]
                    for k in range(0,3):
                        s += '%10.4f ' % origo[j,k]
                    s += "\n"

                fh_dumpL.write(s)

            if(args.dumpG==True):
                mat = t.lcs2mat_4d(lcs,origo,args.cutoff)
                s = stringify(files[ii],jj,mat,sequence[jj],hread=args.read)
                fh_dumpG.write(s)

            if(args.dumpR==True):
                mat = t.lcs2mat_3d(lcs,origo,args.cutoff)
                s = stringify(files[ii],jj,mat,sequence[jj],hread=args.read)
                fh_dumpR.write(s)
 

        
    if(args.dumpG==True):
        fh_dumpG.close()
    if(args.dumpR==True):
        fh_dumpR.close()
    if(args.dumpL==True):
        fh_dumpL.close()


    return 0


