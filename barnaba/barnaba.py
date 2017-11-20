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
