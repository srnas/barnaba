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

#import reader as reader
import mdtraj as md
import btools as bt
import definitions

def stringify(fname,data,seq1,seq2,hread=True):

    if(hread==False):
        data = data.reshape(-1)
        s = "%10s " % (fname)
        s += "".join([("%10.6f " % el) for el in data])
        return s + "\n"
    else:
        s = '# %s \n' % (fname)
        for j in range(data.shape[0]):
            for k in range(data.shape[1]):
                if(all(data[j,k]) == 0.0):
                    continue
                s += '%10s %10s ' % (seq1[j],seq2[k])
                s += "".join(['%10.6f ' % data[j,k,l] for l in range(len(data[j,k]))])
                s+= "\n"
        return s
    
def dump(args):
    

    print "# Calculating DUMP ... (ahah)"

    if(args.dumpG==True):
        fh_dumpG = open(args.name+ ".gvec",'w')

    if(args.dumpR==True):
        fh_dumpR = open(args.name+ ".rvec",'w')
        
    if(args.dumpP==True):
        fh_dumpP = open(args.name+ ".pvec",'w')


    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]

    for i in xrange(0,len(files)):
        
        if(args.pdbs!=None):
            top= files[i]
            cur_pdb = md.load_pdb(files[i])
        else:
            top=args.top
            cur_pdb = md.load_frame(files[i],0,top=args.top)
    
        cur_idx = bt.get_lcs_idx(cur_pdb.topology)


        residues = [cur_pdb.topology.atom(at).residue  for at in cur_idx[0]]
        seq_id = [definitions.residue_dict[rr.name] + str(rr.resSeq)  for rr in residues]
        label = files[i]

        if(len(cur_idx[0])<2):
            print "# Warning: File %s has %d RNA nucleotide(s). Skipping." % (files[i],len(cur_idx[0]))
            continue

        if(args.dumpP == True):
            other_idx = bt.get_other_idx(cur_pdb.topology,args.atom)
            other_residues = [cur_pdb.topology.atom(at).residue  for at in other_idx]
            other_seq_id = [definitions.residue_dict[rr.name] + str(rr.resSeq)  for rr in other_residues]


        # duplicate code - it would be horribly slow otherwise
        if(args.pdbs!=None):
            
            c1 = cur_pdb.xyz[0,cur_idx[0]]
            c2 = cur_pdb.xyz[0,cur_idx[1]]
            c3 = cur_pdb.xyz[0,cur_idx[2]]
            label = files[i]
                    
            if(args.dumpG==True):
                
                mat = bt.get_gmat(c1,c2,c3,args.cutoff)
                s = stringify(label,mat,seq_id,seq_id,hread=args.read)
                fh_dumpG.write(s)
                    
            if(args.dumpR==True):
                mat = bt.get_3dmat_square(c1,c2,c3,args.cutoff)
                s = stringify(label,mat,seq_id,seq_id,hread=args.read)
                fh_dumpR.write(s)
                    
            if(args.dumpP==True):
                c5 = cur_pdb.xyz[0,other_idx] 
                mat = bt.get_other_mat(c1,c2,c3,c5,args.cutoff)
                s = stringify(label,mat,seq_id,other_seq_id,hread=args.read)
                fh_dumpP.write(s)
                
        else:

            # analyze trajectory in chunks of 100
            for chunk in md.iterload(files[i], chunk=100,top=top,stride=args.stride):
                
                for j in range(len(chunk)):
                    
                    c1 = chunk.xyz[j,cur_idx[0]]
                    c2 = chunk.xyz[j,cur_idx[1]]
                    c3 = chunk.xyz[j,cur_idx[2]]
                    
                    label = str(chunk.time[j])
                    
                    if(args.dumpG==True):
                        
                        mat = bt.get_gmat(c1,c2,c3,args.cutoff)
                        s = stringify(label,mat,seq_id,seq_id,hread=args.read)
                        fh_dumpG.write(s)
                        
                    if(args.dumpR==True):
                        mat = bt.get_3dmat_square(c1,c2,c3,args.cutoff)
                        s = stringify(label,mat,seq_id,seq_id,hread=args.read)
                        fh_dumpR.write(s)
                    
                    if(args.dumpP==True):
                        c5 = chunk.xyz[j,other_idx] 
                        mat = bt.get_other_mat(c1,c2,c3,c5,args.cutoff)
                        s = stringify(label,mat,seq_id,other_seq_id,hread=args.read)
                        fh_dumpP.write(s)
                
                                
    if(args.dumpG==True):
        fh_dumpG.close()
    if(args.dumpR==True):
        fh_dumpR.close()
    if(args.dumpP==True):
        fh_dumpP.close()


    return 0


