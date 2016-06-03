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

import mdtraj as md
import btools as bt
import definitions as definitions
import os
import string
import random
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

from scipy.spatial import distance


def snippet(args):

    
    print "# Annotating RNA structures..."

    # check query sequence
    for item in args.seq:
        if(item not in definitions.known_abbrev):
            print "# FATAL Error. Symbol ", item, " not known. Use ACGU NYR"
            return 1
        if(item == "%"):
            print "# Fatal error. Single strand only"
            return 1
        
    fh = open(args.name,'a')

    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]


    for i in xrange(0,len(files)):
        print "#",files[i]

        ii = 0
        name_pref = files[i][0:-4].split("/")[-1]
        
        if(args.pdbs!=None):
            top= files[i]
            # remove aminoacids and heteroatoms from list
            ss = id_generator() + ".pdb"
            cmd = "grep -v \"" + " \| ".join(definitions.aa)+ "\| HETATM \" " +top+" > " + ss
            os.system(cmd)
            try:
                cur_pdb = md.load_pdb(ss)
            except:
                print "# ERROR: could not load ",top, "... skipping"
            os.system("rm " + ss)
        else:
            top=args.top
            cur_pdb = md.load_frame(files[i],0,top=args.top)

        cur_idx = bt.get_lcs_idx(cur_pdb.topology)

        residues = [cur_pdb.topology.atom(at).residue  for at in cur_idx[0]]
        seq_id = [definitions.residue_dict[rr.name] + str(rr.resSeq) + "_" + str(rr.chain.index)  for rr in residues]
        seq = [definitions.residue_dict[rr.name]  for rr in residues]
        res_idxs_full = [rr.index  for rr in residues]
        res_idxs = bt.get_idx(seq,args.seq,bulges=0)
        
        if(len(residues)<len(args.seq)): continue
        if(len(res_idxs)==0): continue
        
        

        # duplicate code - it would be horribly slow otherwise
        if(args.pdbs!=None):
            for el in res_idxs:

                seqs = "".join([seq[p] for p in el ])
                seqs_id = " ".join([seq_id[p] for p in el ])
                resi_full = [res_idxs_full[y] for y in el]
                # skip non-consecutive residues
                if(resi_full[-1]!=resi_full[0]+len(resi_full)-1):
                    continue
                
                tmp_atoms = [atom.index for atom in cur_pdb.topology.atoms if (atom.residue.index in resi_full)]
                trj_slice = cur_pdb[0].atom_slice(tmp_atoms)
                trj_slice.center_coordinates()

                name = name_pref + "_" + seqs + "_" + str(ii).zfill(4) + ".pdb"
                trj_slice.save(name)
                fh.write("%30s %10.4f %d %s \n" % (files[i],0,ii,seqs_id) )
                ii += 1

        else:
            # analyze trajectory in chunks of 100            
            for chunk in md.iterload(files[i], chunk=100,top=top):
                for j in range(len(chunk)):
                    
                    for el in res_idxs:
                        seqs = "".join([seq[p] for p in el ])
                        seqs_id = " ".join([seq_id[p] for p in el ])
                        resi_full = [res_idxs_full[y] for y in el]
                        # skip non-consecutive residues
                        if(resi_full[-1]!=resi_full[0]+len(resi_full)-1):
                            continue
                        tmp_atoms = [atom.index for atom in cur_pdb.topology.atoms if (atom.residue.index in resi_full)]
                        trj_slice = chunk[i].atom_slice(tmp_atoms)
                        trj_slice.center_coordinates()

                        name = name_pref + "_" + seqs + "_" + str(ii).zfill(4) + ".pdb"
                        trj_slice.save(name)
                        fh.write("%30s %10.4f %d %s \n" % (files[i],chunk.time[j],ii,seqs_id) )
                        ii += 1

    fh.close()
    return 0


