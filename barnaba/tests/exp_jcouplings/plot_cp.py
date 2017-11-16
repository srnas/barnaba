import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# sugar

    
def cp_mooren(x):
    return 8*np.cos(x)*np.cos(x)-3.4*np.cos(x) + 0.5

def cp_marino(x):
    return 6.9*np.cos(x)*np.cos(x)-3.4*np.cos(x) + 0.7


fmts = ['o','s','*','d']
cols = ['blue','red','black','yellow']
mooren_err = []
marino_err = []

xx = np.arange(0,2*np.pi,0.01)
ff = np.pi/180.

plt.plot(xx,cp_mooren(xx),c="g",label="mooren")
plt.plot(xx,cp_marino(xx),c="r",label="marino")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("c4p3.txt") if ("#" not in line)])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])
plt.xlabel("Epsilon (rad)")
plt.ylabel("3J (Hz)")
plt.savefig("c4p3.png")
plt.close()
marino_err.extend([cp_marino(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
mooren_err.extend([cp_mooren(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])

plt.plot(xx,cp_mooren(xx),c="g",label="mooren")
plt.plot(xx,cp_marino(xx),c="r",label="marino")
plt.legend()
data = np.array([[float(x) for x in line.split()] for line in open("c4p5.txt") if ("#" not in line)])
data[np.where(data[:,3]<0.0),3] += 360
plt.errorbar(data[:,3]*ff,data[:,1],fmt ='o',xerr=data[:,4]*ff,yerr=data[:,2])
plt.xlabel("Beta (rad)")
plt.ylabel("3J (Hz)")
plt.savefig("c4p5.png")
plt.close()
marino_err.extend([cp_marino(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])
mooren_err.extend([cp_mooren(data[ii,3]*ff)-data[ii,1] for ii in range(data.shape[0])])


print "MARINO",np.sqrt(np.average(np.array(marino_err)**2))
print "MOOREN",np.sqrt(np.average(np.array(mooren_err)**2))



