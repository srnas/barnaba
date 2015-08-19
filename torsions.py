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
import numpy as np
import definitions


def torsions(args):


    print "# Calculating torsion angles..."
    files = args.files

    header = ["# " + str(k) + " " + str(args.__dict__[k]) + "\n" for k in sorted(args.__dict__)]

    if(args.bb):
        fh_bb = open(args.name + ".backbone",'w')
        fh_bb.write("".join(header))
    if(args.pucker):
        fh_pu = open(args.name + ".pucker",'w')
        fh_pu.write("".join(header))
    if(args.jcoupling):
        fh_j3 = open(args.name + ".j3",'w')
        fh_j3.write("".join(header))
        
    for i in xrange(0,len(files)):

        cur_pdb = reader.Pdb(files[i],res_mode=args.res_mode)
        ll = len(cur_pdb.model.sequence)
        if(args.xtc!=None):
            cur_pdb.set_xtc(args.xtc)

        if(args.bb):
            cur_pdb.model.set_bb_index()
        if(args.pucker):
            cur_pdb.model.set_pucker_index()
        if(args.jcoupling):
            cur_pdb.model.set_j3_index()
            
        idx = 0
        eof = True
        while(eof):

            ################################
            ### calculate backbone angles ##
            ################################

            if(args.bb):
                bb_angles, chi_angles = cur_pdb.model.calc_bb_torsion()
                idx1 = 0
                idx2 = 0
                angles = []
                for res in range(ll):
                    string = ""
                    for bbs in range(6):
                        if(cur_pdb.model.missing[res*6+bbs]):
                            string +=  "%10.3f " % float('nan')
                        else:
                            string +=  "%10.3f "%  bb_angles[idx1]
                            idx1 += 1
                            
                    if(cur_pdb.model.chi_missing[res]):
                        string +=  "%10.3f " % float('nan')
                    else:
                        string +=  "%10.3f " %chi_angles[idx2]
                        idx2 += 1
                    angles.append(string)
                    
                # now print to file
                if(args.hread):
                    fh_bb.write("# " + files[i] + " " + str(idx) + "\n")
                    fh_bb.write("#%10s%10s %10s %10s %10s %10s %10s %10s \n" % ("","ALPHA","BETA","GAMMA","DELTA","EPSILON","ZETA","CHI"))
                    for res in range(ll):
                        stri = "%10s " % cur_pdb.model.sequence_id[res]
                        stri += angles[res]
                        fh_bb.write(stri + "\n")
                else:
                    fh_bb.write("%10d %s \n"% (idx,"".join(angles)))
                    

            ################################
            ### calculate pucker angles  # #
            ################################
            if(args.pucker):
                pucker_angles = cur_pdb.model.calc_pucker()
                idx1 = 0
                angles = []
                for res in range(ll):
                    string = ""
                    if(cur_pdb.model.pucker_missing[res]):
                        angles.append("%10.3f " % float('nan'))
                    else:
                        angles.append("".join(["%10.3f "% el for el in pucker_angles[idx1]]))
                        idx1 += 1
                        
                if(args.hread):
                    fh_pu.write("# " + files[i] + " " + str(idx) + "\n")
                    fh_pu.write("#%10s%10s %10s %10s %10s %10s %10s %10s \n" % ("","NU1","NU2","NU3","NU4","NU5","AMPLITUDE","PHASE"))
                    for res in range(ll):
                        stri = "%10s " % cur_pdb.model.sequence_id[res]
                        stri += angles[res]
                        fh_pu.write(stri + "\n")
                else:
                    fh_pu.write("%10d %s \n"% (idx,"".join(angles)))
                            
            ################################
            ### calculate pucker angles  # #
            ################################
            if(args.jcoupling):
                j3_angles = cur_pdb.model.calc_j3()
                idx1 = 0
                angles = []
                lj = len(definitions.j3)
                for res in range(ll):
                    string = ""
                    for bbs in range(lj):
                        if(cur_pdb.model.j3_missing[res*lj+bbs]):
                            string +=  "%10.3f " % float('nan')
                        else:
                            if(args.raw):
                                val = j3_angles[idx1]
                            else:
                                cos = np.cos(j3_angles[idx1])
                                val =  cos*cos*definitions.j3[bbs][2][0] + \
                                       cos*definitions.j3[bbs][2][1] +\
                                       definitions.j3[bbs][2][2]
                            string +=  "%10.3f "%  val
                            idx1 += 1
                    angles.append(string)
                    
                # now print to file
                if(args.hread):
                    fh_j3.write("# " + files[i] + " " + str(idx) + "\n")
                    ss = "#%10s" % ""
                    for el in definitions.j3:
                        ss+= "%10s " % el[0]
                    fh_j3.write(ss + "\n")
                    for res in range(ll):
                        stri = "%10s " % cur_pdb.model.sequence_id[res]
                        stri += angles[res]
                        fh_j3.write(stri + "\n")
                else:
                    fh_j3.write("%10d %s \n"% (idx,"".join(angles)))
                
                        
            
            idx += 1
            if(args.xtc==None):
                eof = cur_pdb.read()
            else:
                eof = cur_pdb.read_xtc()

    if(args.bb):
        fh_bb.close()
    if(args.pucker):
        fh_pu.close()
    if(args.jcoupling):
        fh_j3.close()
    return 0
            
