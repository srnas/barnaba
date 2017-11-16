import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# sugar
def sugar_altona_h1h2(x):
    v = [6.96462,-0.91,1.02629,1.27009,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def sugar_altona_h2h3(x):
    v = [8.28926,-0.91, 0.66772, 0.00193297,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

def sugar_altona_h3h4(x):
    v = [7.96446,-0.91, 0.77241, -0.262475,0]
    cos = np.cos(x+v[4])
    sin = np.sin(x+v[4])
    return v[0]*cos*cos + v[1]*cos + v[2] + v[3]*cos*sin

    
def sugar_davies(x):
    return 10.2*np.cos(x)*np.cos(x)-0.8*np.cos(x)

def sugar_condon(x):
    return 9.67*np.cos(x)*np.cos(x)-2.03*np.cos(x)


fmts = ['o','s','*','d']
cols = ['blue','red','black','yellow']
altona_err = []
davies_err = []
condon_err = []

xx = np.arange(0,2*np.pi,0.01)
ff = np.pi/180.
plt.plot(xx,sugar_altona_h1h2(xx),c="b",label="altona")
plt.plot(xx,sugar_davies(xx),c="g",label="davies")
plt.plot(xx,sugar_condon(xx),c="r",label="condon")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("h1h2.txt")])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])

condon_err.extend([sugar_condon(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
davies_err.extend([sugar_davies(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
altona_err.extend([sugar_altona_h1h2(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])


plt.savefig("h1h2.png")
plt.close()

plt.plot(xx,sugar_altona_h2h3(xx),c="b",label="altona")
plt.plot(xx,sugar_davies(xx),c="g",label="davies")
plt.plot(xx,sugar_condon(xx),c="r",label="condon")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("h2h3.txt")])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])
condon_err.extend([sugar_condon(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
davies_err.extend([sugar_davies(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
altona_err.extend([sugar_altona_h2h3(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])

#print np.sqrt(np.average(np.array(condon_err)**2))
#print np.sqrt(np.average(np.array(davies_err)**2))
#print np.sqrt(np.average(np.array(altona_err)**2))

plt.savefig("h2h3.png")
plt.close()

plt.plot(xx,sugar_altona_h3h4(xx),c="b",label="altona")
plt.plot(xx,sugar_davies(xx),c="g",label="davies")
plt.plot(xx,sugar_condon(xx),c="r",label="condon")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("h3h4.txt")])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])
condon_err.extend([sugar_condon(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
davies_err.extend([sugar_davies(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
altona_err.extend([sugar_altona_h3h4(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])

print "CONDON",np.sqrt(np.average(np.array(condon_err)**2))
print "DAVIES",np.sqrt(np.average(np.array(davies_err)**2))
print "ALTONA",np.sqrt(np.average(np.array(altona_err)**2))


plt.savefig("h3h4.png")
plt.close()

# ERRORS


##################  BETA ############
