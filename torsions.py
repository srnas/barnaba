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
import numpy as np
import mdtraj as md
import btools as bt
import definitions
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import colors,lines


def plot():
    
    # this is for plotting!
    lab1 = [r'$\alpha$',r'$\beta$',r'$\gamma$',r'$\delta$',r'$\epsilon$',r'$\zeta$',r'$\chi$']

    # define colormaps
    cmap = colors.ListedColormap(['white', 'orange', 'green'])
    bounds=[-10,1,20,1000]
    norm = colors.BoundaryNorm(bounds, cmap.N)


    fig = plt.figure()

    # labels for Gauche/Trans
    fig.text(0.17,0.79,"G+")
    line1 = lines.Line2D([np.pi/6,np.pi/6],[7.5,8], lw=0.5, color='k',zorder=10)
    line1.set_clip_on(False)
    line2 = lines.Line2D([np.pi/2,np.pi/2],[7.5,8], lw=0.5, color='k',zorder=10)
    line2.set_clip_on(False)
    
    fig.text(0.295,0.79,"T")
    line3 = lines.Line2D([5*np.pi/6,5*np.pi/6],[7.5,8], lw=0.5, color='k',zorder=10)
    line3.set_clip_on(False)
    line4 = lines.Line2D([7*np.pi/6,7*np.pi/6],[7.5,8], lw=0.5, color='k',zorder=10)
    line4.set_clip_on(False)
    
    fig.text(0.41,0.79,"G-")
    line5 = lines.Line2D([9*np.pi/6,9*np.pi/6],[7.5,8], lw=0.5, color='k',zorder=10)
    line5.set_clip_on(False)
    line6 = lines.Line2D([11*np.pi/6,11*np.pi/6],[7.5,8], lw=0.5, color='k',zorder=10)
    line6.set_clip_on(False)
        

    # lines for A-form helix
    lines_a = [lines.Line2D([-1.08365+2*np.pi,-1.08365+2*np.pi],[0,1], lw=1.5, color='r'),\
               lines.Line2D([-3.13945+2*np.pi,-3.13945+2*np.pi],[0,1], lw=1.5, color='r'),\
               lines.Line2D([0.82790,0.82790],[0,1], lw=1.5, color='r'),\
               lines.Line2D([1.45658,1.45658],[0,1], lw=1.5, color='r'),\
               lines.Line2D([-2.64786+2*np.pi,-2.64786+2*np.pi],[0,1], lw=1.5, color='r'),\
               lines.Line2D([-1.28431+2*np.pi,-1.28431+2*np.pi],[0,1], lw=1.5, color='r'),\
               lines.Line2D([-2.89184+2*np.pi,-2.89184+2*np.pi],[0,1], lw=1.5, color='r')]

    # print histogram a-z
    for pp in range(6):
        ap = np.array(all_angles_b[pp])
        h1, e1 = np.histogram(ap,bins=np.arange(0,2*np.pi,0.1),normed=False)
        ax = plt.subplot2grid((6,2),(pp, 0))
        ax.imshow([h1],extent=(0,2*np.pi,0,1),cmap=cmap, norm=norm,alpha=0.5 )
        
        ax.set_yticks([])
        ax.set_ylabel(lab1[pp],fontsize=14,rotation=0)
        ax.add_line(lines_a[pp])
        if(pp!=5):
            ax.set_xticks([])
        else:
            ax.add_line(line1)
            ax.add_line(line2)
            ax.add_line(line3)
            ax.add_line(line4)
            ax.add_line(line5)
            ax.add_line(line6)

    # print histogram for chi
    ap = np.array(all_angles_b[6])
    h1, e1 = np.histogram(ap,bins=np.arange(0,2*np.pi,0.1),normed=False)
    ax = plt.subplot2grid((6,2),(0, 1))
    ax.imshow([h1],extent=(0,2*np.pi,0,1),cmap=cmap, norm=norm,alpha=0.5 )
    ax.set_yticks([])
    ax.add_line(lines_a[6])
    ax.set_ylabel(lab1[6],fontsize=14,rotation=0)
        

    # polar plot
    ax = plt.subplot2grid((6,2),(3, 1),rowspan=2,aspect='equal')
    z1 = all_angles_s[5]*np.cos(all_angles_s[6])
    z2 = all_angles_s[5]*np.sin(all_angles_s[6])
    hh, ee1,ee2 = np.histogram2d(z1,z2,bins=[np.arange(-1,1,0.05),np.arange(-1,1,0.05)],normed=False)

    ax.pcolor(ee1,ee2,hh.T,cmap=cmap, norm=norm)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.axis('off')
    circle2=plt.Circle((0,0),0.99,edgecolor="k",fill=False)
    ax.add_artist(circle2)
    
    dd = 1.17
    dx = 0.15
    dy = 0.06
    labs1 = ["0$^\circ$","36$^\circ$","72$^\circ$","108$^\circ$","144$^\circ$"]
    labs2 = ["180$^\circ$","216$^\circ$","252$^\circ$","288$^\circ$","324$^\circ$"]
    for jj in range(5):
        xx = np.cos(jj*np.pi/5)
        yy = np.sin(jj*np.pi/5)
        line_a1 = lines.Line2D([-xx,xx],[-yy,yy], lw=0.5, color='k',ls='dotted')
        ax.add_line(line_a1)
        ax.text(xx*dd-dx,yy*dd-dy,labs1[jj])
        ax.text(-xx*dd-dx,-yy*dd-dy,labs2[jj])
        
    fig.subplots_adjust(hspace=-0.7)
    plt.savefig("pippo.pdf")


def stringify(angles,miss,seq,hread):

    dd = angles.shape[1]/len(seq)
    result = []
    for t in xrange(angles.shape[0]):
        string = ""
        for ii in xrange(len(seq)):

            i1 = dd*ii
            stri = "".join(["%10.6f "  % angles[t,i1+jj] if ([jj,ii] not in miss) else "%10s " % ("NaN") for jj in range(dd) ]) 

            # calculate pucker angles
            x1 = angles[t,i1+11] + angles[t,i1+8] - angles[t,i1+10] - angles[t,i1+7]
            x2 = 3.0776835*angles[t,i1+9]
            p0 = np.arctan2(x1,x2)
            if(p0<0.0): p0+=2.0*np.pi
            tm = angles[t,i1+9]/np.cos(p0)

            stri += "%10.6f %10.6f" % (p0,tm)
            if(hread):
                stri = "%10s %s \n " % (seq[ii], stri)
            string += stri
        result.append(string)
        
    return result

def dihedral_idx(topology):

    # create and array with mising atoms
    miss = []
    indeces = []
    seq = []

    for ii,rr in enumerate(topology.residues):
        if(rr.name in definitions.rna):

            idxs = [[0,0,0,0] for x in range(12)]
            seq.append("%s_%d_%d" % (rr.name,rr.resSeq,rr.chain.index))
            jj = len(seq)-1
            
            # alpha
            if(ii!=0):
                rr_m = topology.residue(ii-1)
                # check chain and atoms
                if(rr_m.index +1 == rr.index and rr_m.chain.index == rr.chain.index):
                    try: idxs[0] = [rr_m.atom("O3'").index, rr.atom("P").index, rr.atom("O5'").index, rr.atom("C5'").index]
                    except:
                        miss.append([0,jj])
                else:
                    miss.append([0,jj])
            else:
                miss.append([0,jj])
                
            # beta
            try: idxs[1] = [rr.atom("P").index, rr.atom("O5'").index, rr.atom("C5'").index, rr.atom("C4'").index]
            except: miss.append([1,jj])
            # gamma
            try: idxs[2] = [rr.atom("O5'").index, rr.atom("C5'").index, rr.atom("C4'").index, rr.atom("C3'").index]
            except: miss.append([2,jj])
            # delta
            try: idxs[3] = [rr.atom("C5'").index, rr.atom("C4'").index, rr.atom("C3'").index, rr.atom("O3'").index]
            except: miss.append([3,jj])

            # epsilon 
            try:
                rr_p = topology.residue(ii+1)
                if(rr_p.index - 1 == rr.index and rr_p.chain.index == rr.chain.index):
                    try: idxs[4]  = [rr.atom("C4'").index, rr.atom("C3'").index, rr.atom("O3'").index, rr_p.atom("P").index]
                    except: miss.append([4,jj])
                else: miss.append([4,jj])
            except: miss.append([4,jj])

                
            # zeta
            try:
                rr_p = topology.residue(ii+1)
                if(rr_p.index - 1 == rr.index and rr_p.chain.index == rr.chain.index):
                    try: idxs[5]  = [rr.atom("C3'").index, rr.atom("O3'").index, rr_p.atom("P").index, rr_p.atom("O5'").index]
                    except: miss.append([5,jj])
                else: miss.append([5,jj])
            except: miss.append([5,jj])
                
                
            if(definitions.pypu_dict[rr.name] == "R"):
                try: idxs[6]  =  [rr.atom("O4'").index, rr.atom("C1'").index, rr.atom("N9").index, rr.atom("C4").index]
                except: miss.append([6,jj])
            if(definitions.pypu_dict[rr.name] == "Y"):
                try: idxs[6]  =  [rr.atom("O4'").index, rr.atom("C1'").index, rr.atom("N1").index, rr.atom("C2").index]
                except: miss.append([6,jj])
                
            # sugar 
            try: idxs[7]  = [rr.atom("C4'").index, rr.atom("O4'").index, rr.atom("C1'").index, rr.atom("C2'").index]
            except: miss.append([7,jj])
            
            try: idxs[8]  = [rr.atom("O4'").index, rr.atom("C1'").index, rr.atom("C2'").index, rr.atom("C3'").index]
            except: miss.append([8,jj])
            
            try: idxs[9]  = [rr.atom("C1'").index, rr.atom("C2'").index, rr.atom("C3'").index, rr.atom("C4'").index]
            except: miss.append([9,jj])
            
            try: idxs[10]  = [rr.atom("C2'").index, rr.atom("C3'").index, rr.atom("C4'").index, rr.atom("O4'").index]
            except: miss.append([10,jj])
            
            try: idxs[11]  = [rr.atom("C3'").index, rr.atom("C4'").index, rr.atom("O4'").index, rr.atom("C1'").index]
            except: miss.append([11,jj])        
            
            indeces.extend(idxs)

    return indeces,seq,miss

def torsions(args):

    print "# Calculating torsion angles..."
    fh = open(args.name,'a')
    labs = ["# ","ALPHA","BETA","GAMMA","DELTA","EPSILON","ZETA","CHI","NU1","NU2","NU3","NU4","NU5","AMPLITUDE","PHASE"]
    fh.write("".join(["%10s " % la for la in labs]) + "\n")
    
    if(args.pdbs!=None):
        files = args.pdbs
    else:
        files = [args.trj]

    
    for i in xrange(0,len(files)):
        print "#",files[i]

        if(args.pdbs!=None):
            cur_pdb = md.load_pdb(files[i])
            idxs,seq,miss = dihedral_idx(cur_pdb.topology)
            angles = md.compute_dihedrals(cur_pdb,idxs,periodic=False)
            string = stringify(angles,miss,seq,args.hread)
            if(args.hread):
                fh.write("# file %s \n" % (files[i]))
                fh.write(string[0])
            else:
                fh.write("%s %s \n" % (files[i],string[0]))
                
        else:
            cur_pdb = md.load_frame(files[i],0,top=args.top)
            idxs,seq,miss = dihedral_idx(cur_pdb.topology)
            
            for chunk in md.iterload(files[i], chunk=100,top=args.top):
                angles = md.compute_dihedrals(cur_pdb,idxs,periodic=False)
                string = stringify(angles,miss,seq,args.hread)
                
                for t in range(len(string)):
                    if(args.hread):
                        fh.write("# file %s time=%10f \n" % (files[i], chunk[t].time))
                        fh.write(string[t])
                    else:
                        fh.write("%10f %s \n" % (chunk[t].time,string[t]))
                
   

        
            
    return 0
            

