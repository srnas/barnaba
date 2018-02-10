from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cp = sns.color_palette(n_colors=9)
sns.set_style("white")
sns.set_context("paper")

# sugar
def sugar_hasnoot_h1h2(x):
    v = [6.96462,-0.91,1.02629,1.27009,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def sugar_hasnoot_h2h3(x):
    v = [8.28926,-0.91, 0.66772, 0.00193297,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def sugar_hasnoot_h3h4(x):
    v = [7.96446,-0.91, 0.77241, -0.262475,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def hasnoot_h4h5p(x):
    v = [8.313,-0.99,1.373,0.27,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def hasnoot_h4h5s(x):
    v = [8.313,-0.99,1.373,-4.752,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin
    
def sugar_davies(x):
    return 10.2*np.cos(x)*np.cos(x)-0.8*np.cos(x)

def sugar_condon(x):
    return 9.67*np.cos(x)*np.cos(x)-2.03*np.cos(x)

def cp_mooren(x):
    return 8*np.cos(x)*np.cos(x)-3.4*np.cos(x) + 0.5

def cp_marino(x):
    return 6.9*np.cos(x)*np.cos(x)-3.4*np.cos(x) + 0.7

def hp_lankhorst(x):
    #x += 
    return 15.3*np.cos(x)*np.cos(x)-6.1*np.cos(x) + 1.6

def hp_mooren(x):
    return 15.3*np.cos(x)*np.cos(x)-6.2*np.cos(x) + 1.5

def hp_lee(x):
    return 18.1*np.cos(x)*np.cos(x)-4.8*np.cos(x) + 1.5
def ippel_c42(x):
    return 4.7*np.cos(x)*np.cos(x)+2.3*np.cos(x) + 0.1

   
def ippel_c86(x):
    return 4.5*np.cos(x)*np.cos(x)-0.6*np.cos(x) + 0.1

def hh_davies(x):
    return 9.7*np.cos(x)*np.cos(x)-1.8*np.cos(x)


xx = np.arange(0,2*np.pi,0.01)
ff = np.pi/180.
errors = []
labels = []
fig,ax = plt.subplots(5,3,figsize=(8.3,8),sharex=True)

dataj = np.array([[float(x) for x in line.split()] for line in open("jcouplings.dat") if len(line.split())>0])
dataa = np.array([[float(x) for x in line.split()] for line in open("jcouplings_raw.dat") if len(line.split())>0])

print "H1H2"
diff = np.abs(sugar_davies(dataa[:,0])-dataj[:,0])
ax[0,0].scatter(np.arange(dataa.shape[0]),diff,s=1.0)
ax[0,0].set_ylim(0,1.2*np.max(diff))
print "%d %.1e %.1e " % (0, np.max(diff),np.median(diff))

print "H2H3"
diff = np.abs(sugar_davies(dataa[:,1])-dataj[:,1])
ax[0,1].scatter(np.arange(dataa.shape[0]),diff,s=1.0)
ax[0,1].set_ylim(0,1.2*np.max(diff))
print "%d %.1e %.1e " % (1, np.max(diff),np.median(diff))

bb=np.pi*(2./3.)

print "H3H4"
diff = np.abs(sugar_davies(dataa[:,2])-dataj[:,2])
ax[0,2].scatter(np.arange(dataa.shape[0]),diff,s=1.0)
ax[0,2].set_ylim(0,1.2*np.max(diff))
print "%d %.1e %.1e " % (2, np.max(diff),np.median(diff))


#ax[0,0].set_ylim(0.0001,1.2*np.max(diff))
#collections.OrderedDict([("H1H2",0),("H2H3",1),("H3H4",2),\
#                             ("1H5P",3),("2H5P",3),("C4Pb",3),\
#                             ("1H5H4",4),("2H5H4",4),\
#                             ("H3P",5),("C4Pe",5),\
#                             ("H1C2/4",6),("H1C6/8",6)])

print "1H5P"
diff = np.abs(hp_lee(dataa[:,3]-bb)-dataj[:,3])
diff = diff[np.isfinite(diff)]
ax[2,0].scatter(np.arange(diff.shape[0]),diff,s=1.0)
print "%d %.1e %.1e " % (3, np.max(diff),np.median(diff))
ax[2,0].set_ylim(0,1.2*np.max(diff))

print "2H5P"
diff = np.abs(hp_lee(dataa[:,3]+bb)-dataj[:,4])
diff = diff[np.isfinite(diff)]
ax[2,1].scatter(np.arange(diff.shape[0]),diff,s=1.0)
print "%d %.1e %.1e " % (4, np.max(diff),np.median(diff))
ax[2,1].set_ylim(0,1.2*np.max(diff))

print "C4Pb"
diff = np.abs(cp_marino(dataa[:,3])-dataj[:,5])
diff = diff[np.isfinite(diff)]
ax[1,0].scatter(np.arange(diff.shape[0]),diff,s=1.0)
ax[1,0].set_ylim(0,1.2*np.max(diff))
print "%d %.1e %.1e " % (2, np.max(diff),np.median(diff))


print "1H5H4"
diff = np.abs(hasnoot_h4h5p(dataa[:,4]-bb)-dataj[:,6])
diff = diff[np.isfinite(diff)]
ax[3,0].scatter(np.arange(diff.shape[0]),diff,s=1.0)
ax[3,0].set_ylim(0,1.2*np.max(diff))
print "%d %.1e %.1e " % (2, np.max(diff),np.median(diff))

print "2H5H4"
diff = np.abs(hasnoot_h4h5s(dataa[:,4])-dataj[:,7])
diff = diff[np.isfinite(diff)]
ax[3,1].scatter(np.arange(diff.shape[0]),diff,s=1.0)
ax[3,1].set_ylim(0,1.2*np.max(diff))
print "%d %.1e %.1e " % (2, np.max(diff),np.median(diff))


print "H3P"
diff = np.abs(hp_lankhorst(dataa[:,5]+bb)-dataj[:,8])
diff = diff[np.isfinite(diff)]
ax[2,2].scatter(np.arange(diff.shape[0]),diff,s=1.0)
print "%d %.1e %.1e " % (4, np.max(diff),np.median(diff))
ax[2,2].set_ylim(0,1.2*np.max(diff))


print "C4Pe"
diff = np.abs(cp_marino(dataa[:,5])-dataj[:,9])
diff = diff[np.isfinite(diff)]
ax[1,1].scatter(np.arange(diff.shape[0]),diff,s=1.0)
print "%d %.1e %.1e " % (4, np.max(diff),np.median(diff))
ax[1,1].set_ylim(0,1.2*np.max(diff))


print "chi1"
diff = np.abs(ippel_c42(dataa[:,6]-0.5*bb)-dataj[:,10])
diff = diff[np.isfinite(diff)]
ax[4,0].scatter(np.arange(diff.shape[0]),diff,s=1.0)
print "%d %.1e %.1e " % (4, np.max(diff),np.median(diff))
ax[4,0].set_ylim(0,1.2*np.max(diff))

print "chi2"
diff = np.abs(ippel_c86(dataa[:,6]-0.5*bb)-dataj[:,11])
diff = diff[np.isfinite(diff)]
ax[4,1].scatter(np.arange(diff.shape[0]),diff,s=1.0)
print "%d %.1e %.1e " % (4, np.max(diff),np.median(diff))
ax[4,1].set_ylim(0,1.2*np.max(diff))



fig.savefig("karplus_chcek.png",dpi=400)


exit()

