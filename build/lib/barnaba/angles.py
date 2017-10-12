import mdtraj as md
import numpy as np
import sys
import nucleic
import definitions

    
    
def backbone_angles(filename,topology=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    rna_seq,angles  = backbone_angles_traj(traj)
    return rna_seq,angles



    
    
def sugar_angles(filename,topology=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    rna_seq,angles  = sugar_angles_traj(traj)
    return rna_seq,angles
