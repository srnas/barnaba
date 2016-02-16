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

import numpy as np
import definitions
import mdtraj as md
import btools as bt
import re

####################### MOTIF #########################



def noe(args):

    print "# Calculating NOE distances..."
    fh = open(args.name,'a')

    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]
        
    for i in xrange(0,len(files)):
        print "#",files[i]

        if(args.pdbs!=None):
            top= files[i]
            cur_pdb = md.load_frame(files[i],0,top=files[i])
        else:
            top=args.top
            cur_pdb = md.load_frame(files[i],0,top=args.top)

        # indexes of all hydrogens
        rna_idxs = bt.get_rna(cur_pdb.topology)

        idxs = [ii for ii in rna_idxs if ((re.match('[1-2]H', cur_pdb.topology.atom(ii).name) is not None) or (re.match('H', cur_pdb.topology.atom(ii).name) is not None))]

        assert len(idxs)>1, "# Fatal error. No Hydrogen atoms in file %s" % files[i]
        
        pairs = []
        pairs_lab = []
        for i1 in range(len(idxs)):
            for i2 in range(i1+1,len(idxs)):
                pairs.append([idxs[i1],idxs[i2]])
                pairs_lab.append([cur_pdb.topology.atom(idxs[i1]),cur_pdb.topology.atom(idxs[i2])])

        if(args.pdbs!=None):
            data = md.compute_distances(cur_pdb,np.array(pairs),periodic=False)
            
        else:
            data = []
            for chunk in md.iterload(files[i], chunk=100,top=top):
                dd = md.compute_distances(chunk,np.array(pairs),periodic=False)
                data.extend(dd)
            data = np.array(data)

        
        data6 = np.power(data,-6)
        avg = np.power(np.average(data6,axis=0),-1./6.)
        idx_m = np.where(avg<args.cutoff)
        
        bins = args.nbins

        blocks = np.linspace(0,len(data),bins)
        blocks = [int(el) for el in blocks]

        string = ""
        for i1 in idx_m[0]:

            block_avg = [np.power(np.average(data6[blocks[k]:blocks[k+1],i1]),-1./6.) for k in range(len(blocks)-1)]
            block_error = np.std(block_avg)/np.sqrt(bins-1)
            string += "%12s %12s %10.6f %10.6f \n" % (pairs_lab[i1][0],pairs_lab[i1][1], avg[i1] ,block_error)
        fh.write(string)

    fh.close()
