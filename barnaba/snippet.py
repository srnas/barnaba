#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2017 Sandro Bottaro (sandro.bottaro@bio.ku.dk)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import definitions as definitions
import reader as  reader
import numpy as np

atoms = ["C2","C4","C6"]
      
        
def snippet(pdb,sequence):
    
   # check query sequence
   for item in sequence:
      if(item not in definitions.known_abbrev):
         print "# FATAL Error. Symbol ", item, " not known. Use ACGU NYR"
         return 1
      if(item == "%"):
         print "# Fatal error. Single strand only"
         return 1

   ll = [len(el) for el in sequence]
   cur_pdb = reader.Pdb(pdb,res_mode="R",permissive=True)
   cur_len = len(cur_pdb.model.sequence)
   indeces = definitions.get_idx(cur_pdb.model.sequence,sequence,bulges=0)
   
         
   idx = 0
   ii = 0
   while(idx>=0):
      
      for index in indeces:

         # do checks
         skip = False
         for k,res in enumerate(index):
            rr =  cur_pdb.model[res]
            # check that atoms in the base are in place
            out = [rr.get_idx(aa) for aa in atoms]
            if(np.isnan(np.sum(out))):
               skip = True

            # check connectivity
            if(k<len(index)-1):
               # check that O3' and P are connected
               rrp = cur_pdb.model[index[k+1]]
               dd = np.sqrt(np.sum(( np.array(rrp["P"]) -np.array(rr["O3'"]))**2))
               if(dd>1.7):
                  skip = True
         if(skip):
            continue
            
         name_pref = pdb[0:-4].split("/")[-1]
         new_pdb = "%s_%s_%05d.pdb" % (name_pref,cur_pdb.model.sequence_id[index[0]],ii)
         
         fh_pdb = open(new_pdb,'w')
         fh_pdb.write(cur_pdb.model.string_pdb(index,noP=True,center=True))
         fh_pdb.close()
         
         ii += 1
      idx = cur_pdb.read()
      
   return 0


