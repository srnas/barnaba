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

import reader as reader

def stringify(fname,jj,data,seq,hread=True):

    if(hread==False):
        data = data.reshape(-1)
        s = ''
        for el in data:
            s += '%10.6f ' % el
        s += '%s.%i \n' % (fname,jj)
        return s
    else:
        s = '# %s %i \n' % (fname,jj)
        for j in range(data.shape[0]):
            for k in range(data.shape[0]):
                #if(all(data[j,k]) != 0.0):
                s += '%10s %10s ' % (seq[j],seq[k])
                for l in range(len(data[j,k])):
                    s += '%10.6f ' % data[j,k,l]
                s+= "\n"
        return s
    
def dump(args):
    

    files = args.files
    print "# Calculating DUMP ... (ahah)"

    if(args.dumpG==True):
        fh_dumpG = open(args.name+ ".gvec",'w')

    if(args.dumpR==True):
        fh_dumpR = open(args.name+ ".rvec",'w')
        
    if(args.dumpP==True):
        fh_dumpP = open(args.name+ ".pvec",'w')


    for i in xrange(0,len(files)):
        
        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)
        cur_len = len(cur_pdb.model.sequence)
        assert cur_len>2, "# Fatal error: less than 2 nucleotides %d \n" %(cur_len)
        if(args.xtc!=None):
            cur_pdb.set_xtc(args.xtc)
            
        idx = 0
        eof = True
        while(eof):
            if(args.dumpG==True):
                mat = cur_pdb.model.get_gmat(args.cutoff)
                s = stringify(files[i],idx,mat,cur_pdb.model.sequence_id,hread=args.read)
                fh_dumpG.write(s)
            if(args.dumpR==True):
                mat = cur_pdb.model.get_3dmat_square(args.cutoff)
                s = stringify(files[i],idx,mat,cur_pdb.model.sequence_id,hread=args.read)
                fh_dumpR.write(s)
            if(args.dumpP==True):
                mat = cur_pdb.model.get_other_mat(args.cutoff,args.atomtype)
                s = stringify(files[i],idx,mat,cur_pdb.model.sequence_id,hread=args.read)
                fh_dumpP.write(s)
                
            idx += 1
            if(args.xtc==None):
                eof = cur_pdb.read()
            else:
                eof = cur_pdb.read_xtc()
                                
    if(args.dumpG==True):
        fh_dumpG.close()
    if(args.dumpR==True):
        fh_dumpR.close()
    if(args.dumpP==True):
        fh_dumpP.close()


    return 0


