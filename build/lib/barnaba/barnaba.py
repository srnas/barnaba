import sys
import functions
import mdtraj as md

def ermsd(reference,target,cutoff=2.4,topology=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
        
    warn += "# Loaded target %s \n" % target
    sys.stderr.write(warn)

    return functions.ermsd_traj(ref,traj)


def dump_rvec(filename,topology=None,cutoff=2.4):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return  functions.dump_rvec_traj(traj,cutoff)


def dump_gvec(filename,topology=None,cutoff=2.4):
    
    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)

    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.dump_gvec_traj(traj,cutoff)


def rmsd(reference,target,topology=None,out=None):

    ref = md.load(reference)
    warn =  "# Loaded reference %s \n" % reference
        
    if(topology==None):
        traj = md.load(target)
    else:
        traj = md.load(target,top=topology)
    warn += "# Loaded target %s \n" % target

    return functions.rmsd_traj(ref,traj,out)


def backbone_angles(filename,topology=None,residue=None,angle=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.backbone_angles_traj(traj,residue,angle)


    
def sugar_angles(filename,topology=None,residues=None):

    if(topology==None):
        traj = md.load(filename)
    else:
        traj = md.load(filename,top=topology)
    warn = "# Loading %s \n" % filename
    sys.stderr.write(warn)
    return functions.sugar_angles_traj(traj,residues)
