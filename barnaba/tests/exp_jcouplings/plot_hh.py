import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
cp = sns.color_palette()
# sugar
def altona_h4h5p(x):
    v = [8.313,-0.99,1.373,0.27,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def altona_h4h5s(x):
    v = [8.313,-0.99,1.373,-4.752,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin
   
def hh_hruska(x):
    return 9.7*np.cos(x)*np.cos(x)-1.8*np.cos(x)



fmts = ['o','s','*','d']
cols = ['blue','red','black','yellow']
altona_err = []
hruska_err = []

xx = np.arange(0,2*np.pi,0.01)
ff = np.pi/180.
plt.plot(xx,altona_h4h5p(xx-2.094395),c="b",label="altona")
plt.plot(xx,hh_hruska(xx-2.094395),c="r",label="condon")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("h4h5p.txt")])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])
#condon_err.extend([sugar_condon(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
#davies_err.extend([sugar_davies(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
altona_err.extend([altona_h4h5p(data[ii,3]*ff-2.094395)-data[ii,1] for ii in range(data.shape[0])])
hruska_err.extend([hh_hruska(data[ii,3]*ff-2.094395)-data[ii,1] for ii in range(data.shape[0])])
plt.savefig("h4h4p.png")
plt.close()

plt.plot(xx,altona_h4h5s(xx),c="b",label="altona")
plt.plot(xx,hh_hruska(xx),c="r",label="condon")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("h4h5s.txt")])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])
#condon_err.extend([sugar_condon(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
#davies_err.extend([sugar_davies(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
altona_err.extend([altona_h4h5s(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
hruska_err.extend([hh_hruska(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
plt.savefig("h4h4s.png")
plt.close()

#print "CONDON",np.sqrt(np.average(np.array(condon_err)**2))
#print "DAVIES",np.sqrt(np.average(np.array(davies_err)**2))
#print "ALTONA",np.sqrt(np.average(np.array(altona_err)**2))



# ERRORS


##################  BETA ############
