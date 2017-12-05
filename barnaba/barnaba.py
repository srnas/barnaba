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

import sys
import functions
import mdtraj as md
import definitions

def ermsd(reference,target,cutoff=2.4,topology=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
        
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return functions.ermsd_traj(ref,traj,cutoff=cutoff)


def dump_rvec(filename,topology=None,cutoff=2.4):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return  functions.dump_rvec_traj(traj,cutoff=cutoff)


def dump_gvec(filename,topology=None,cutoff=2.4):
    
    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.dump_gvec_traj(traj,cutoff=cutoff)

def annotate(filename,topology=None):
    
    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.annotate_traj(traj)


def rmsd(reference,target,topology=None,out=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target

    return functions.rmsd_traj(ref,traj,out=out)


def backbone_angles(filename,topology=None,residues=None,angles=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.backbone_angles_traj(traj,residues=residues,angles=angles)


    
def sugar_angles(filename,topology=None,residues=None,angles=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.sugar_angles_traj(traj,residues=residues,angles=angles)

def pucker_angles(filename,topology=None,residues=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.pucker_angles_traj(traj,residues=residues)


def jcouplings(filename,topology=None,residues=None,couplings=None,raw=False):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.jcouplings_traj(traj,residues=residues,couplings=couplings,raw=raw)


def ss_motif(query,target,topology=None,treshold=0.8,cutoff=2.4,sequence=None,out=None,bulges=0):

    ref = md.load(query)
    warn =  "# Loaded query %s \n" % query
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)
    
    return functions.ss_motif_traj(ref,traj,treshold=treshold,cutoff=cutoff,sequence=sequence,out=out,bulges=bulges)


def ds_motif(query,target,l1,l2,treshold=0.9,cutoff=2.4,topology=None,sequence=None,bulges=0,out=None):

    ref = md.load(query)
    warn =  "# Loaded query %s \n" % query        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return functions.ds_motif_traj(ref,traj,l1,l2,treshold=treshold,sequence=sequence,bulges=bulges,out=out)


## DOT-BRACKET ##
def dot_bracket(pairings,sequence):

    ll = len(sequence)
    dot_bracket = []
    for k,pp in enumerate(pairings):
        openings = []
        closings = []
        
        for e in range(len(pp[0])):
            if(pp[1][e]=="WCc"):
                idx1 = pp[0][e][0]
                idx2 = pp[0][e][1]
                if(idx1 in openings):
                    warn = "# Frame %d, residue %s has a double WC base-pairs. " % (k,sequence[idx1])
                    warn += " Dot-bracket annotation set to xxx \n"
                    sys.stderr.write(warn)
                    dot_bracket.append("".join(["x"]*ll))                            
                    continue
                
                if(idx2 in closings):
                    warn = "# Frame %d, residue %s has a double WC base-pairs. " % (k,sequence[idx2])
                    warn += " Dot-bracket annotation set to xxx \n"
                    sys.stderr.write(warn)
                    dot_bracket.append("".join(["x"]*ll))
                    continue
                
                openings.append(idx1)
                closings.append(idx2)

        # check pseudoknots
        dotbr = ['.']*ll
        levels = [-1]*ll
        for idx1 in xrange(len(openings)):
            start1 = openings[idx1]
            end1 = closings[idx1]
            # up one level
            levels[start1] += 1
            levels[end1] += 1
        
            for idx2 in xrange(len(closings)):
                end2 = closings[idx2]
                if(levels[end2] == levels[start1]):
                    if(end2 > start1 and end2 < end1):
                        levels[start1] += 1
                        levels[end1] += 1
                        
        for idx1 in xrange(len(openings)):
            start1 = openings[idx1]
            end1 = closings[idx1]
            
            dotbr[start1] = definitions.op[levels[start1]]
            dotbr[end1] = definitions.cl[levels[end1]]
        dot_bracket.append("".join(dotbr))
    return dot_bracket


