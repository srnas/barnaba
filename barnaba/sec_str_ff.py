#!/usr/bin/env python

import re, sys, math

import numpy as np
from numpy.linalg import *
from numpy.random import *
from sec_str_constants import *
import  datetime  


def energy(pos, param, write_energy=False):
    E_i = 0
    for ind, par in param.items():
        if ind == 0:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k, d = par[:,3], par[:,4]
            n_p12 = norm(p1-p2, axis=1)
            E_i += .5 * np.sum(k * (n_p12 - d)**2)
        elif ind == 1:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k, d = par[:,3], par[:,4]
            n_p12 = norm(p1-p2, axis=1)
            E_i += .5 * np.sum(np.heaviside(d-n_p12, 1) * k * (n_p12 - d)**2)
        elif ind == 2:
            i1, i2, i3 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int)
            p1, p2, p3 = pos[i1], pos[i2], pos[i3]
            k, a0 = par[:,4], par[:,5]
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            E_i += 0.5 * np.sum(k * (angle - a0)**2)
        elif ind == 4:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k = par[:,3]
            _p3 = [p2[:,0]+30, 0.5*(p1[:,1]+p2[:,1])]
            p3 = np.swapaxes(_p3,0,1)
            a0 = 0
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            E_i += 0.5 * np.sum(k * (angle - a0)**2)
        elif ind == 5:
            i1, i2, i3 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int)
            p1, p2, p3 = pos[i1], pos[i2], pos[i3]
            k, a0 = par[:,4], par[:,5]
            a02 = np.full((k.size), np.pi*.5)
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            E_i += 0.5 * np.sum(np.heaviside(a02-angle, 1) * k * (angle - a0)**2)
        elif ind == 6:
            i1, i2, i3, i4 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int), par[:,4].astype(int)
            p1, p2, p3, p4 = pos[i1], pos[i2], pos[i3], pos[i4]
            k, a0 = par[:,5], par[:,6]
            n_p12, n_p34 = norm(p1-p2, axis=1), norm(p4-p3, axis=1)
            arg = np.multiply(p2-p1, p4-p3).sum(1)/n_p12/n_p34
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            E_i += 0.5 * np.sum(k * (angle - a0)**2)
        elif ind == 7:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k = par[:,3]
            n_p12 = norm(p1-p2, axis=1)
            E_i += np.sum(k / n_p12)

    return E_i

#def force(pos, param, write_force):
def force(pos, param, __i1, __i2, __k_rep, __d_rep, __k_rep_lr, write_force):
    _shape = pos.shape
    n = _shape[0]
    F = np.zeros(_shape)
    for ind, par in param.items():
        if ind == 0:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            p12 = p1-p2
            k, d = np.column_stack((par[:,3], par[:,3])), np.column_stack((par[:,4], par[:,4]))
            n_p12 = norm(p12, axis=1)
            n12 = np.column_stack((n_p12, n_p12))
            diff = n12 - d
            f = np.zeros((n, n, 2))
            f[i1, i2] = -  k * diff * p12 * np.power(n12, -1)
            F += np.sum(f, axis=1)
            F -= np.sum(f, axis=0)

  #          f1 = f
  #          f2 = -f
  #          if np.unique(i1).size == i1.size and np.unique(i2).size == i2.size:
  #              F[i1] += f1
  #              F[i2] += f2
  #          else:    
  #              s = np.unique(np.append(i1,i2))
  #              for i in s:
  #                  F[i] += np.sum(np.append(f1[np.where(i1==i)], f2[np.where(i2==i)], axis=0), axis=0) 
  #             #     print("f1+f2", i, np.sum(np.append(f1[np.where(i1==i)], f2[np.where(i2==i)], axis=0), axis=0))
  #            #      print("     f2", np.sum(f2[np.where(i2==i)], axis=0))
#      #      print(f)        
#            f2 = np.zeros((n, n, 2))
#            f2[i1, i2] = - k * diff * p12 * np.power(n12, -1)
#      #      print(f2[i1, i2])
#            print(np.sum(f2, axis=1)-np.sum(f2, axis=0))
#         #   print(-np.sum(f2, axis=0))

        elif ind == 1:
            i1, i2 = __i1, __i2
            p1, p2 = np.take(pos, i1, axis=0), np.take(pos, i2, axis=0)
            p12 = p1-p2
            k = __k_rep
            d = __d_rep
            n_p12 = norm(p12, axis=1)
            diff = n_p12 - d[:,0]
            f = np.zeros((n, n, 2))
            _i  = np.where(diff < 0)
            i1, i2 = i1[_i], i2[_i]
            p12 = p12[_i]
            k = k[_i]
            n12 = np.column_stack((n_p12[_i], n_p12[_i]))
            diff = n12 - d[_i]


            f[i1, i2] = - k * diff * p12 * np.power(n12, -1)
        #    print("f 16,17", f[16, 17], norm(pos[16]-pos[17]))
            F += np.sum(f, axis=1)
            F -= np.sum(f, axis=0)
        elif ind == 2:
        #    print(par)
            i1, i2, i3 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int)
            ii = np.array([i1, i2, i3])
            p1, p2, p3 = pos[i1], pos[i2], pos[i3]
            k, a0 = par[:,4], par[:,5]
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-7 
            arg[np.where(arg <= -1.)] = -1.+1e-7 
            angle = np.arccos(arg)
        #    print(angle)
            _f = - k * (angle - a0) * (-1./np.sqrt(1.-(arg)**2))
            for comp in range(2):
                f1 = _f * ( (p2[:,comp]-p3[:,comp])/n_p12/n_p23 + (p2[:,comp]-p1[:,comp])*arg/n_p12**2 )
                f2 = _f * ( (p1[:,comp]-2*p2[:,comp]+p3[:,comp])/n_p12/n_p23 + (p3[:,comp]-p2[:,comp])*arg/n_p23**2 -
                        (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
                f3 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p23 - (p3[:,comp]-p2[:,comp])*arg/n_p23**2 )
                ff = np.array([f1, f2, f3])
                for j in range(3):
                    if np.unique(ii[j]).size < ii[j].size:
                        for ik, _i in np.ndenumerate(ii[j]):
                            F[_i][comp] += ff[j][ik]
                    else:
                        F[ii[j], comp] += ff[j]
        elif ind == 4:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            ii = [i1, i2]
            p1, p2 = pos[i1], pos[i2]
            k = par[:,3]
            _p3 = [p2[:,0]+30, 0.5*(p1[:,1]+p2[:,1])]
            p3 = np.swapaxes(_p3,0,1)
            a0 = np.zeros((k.size))
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            _f = - k * (angle - a0) * (-1./np.sqrt(1.-(arg)**2))
            for comp in range(2):
                f1 = _f * ( (p2[:,comp]-p3[:,comp])/n_p12/n_p23 + (p2[:,comp]-p1[:,comp])*arg/n_p12**2 )
                f2 = _f * ( (p1[:,comp]-2*p2[:,comp]+p3[:,comp])/n_p12/n_p23 + (p3[:,comp]-p2[:,comp])*arg/n_p23**2 -
                        (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
                ff = [f1, f2]
                for j in range(2):
                    if len(set(ii[j])) < ii[j].size:
                        for ik, _i in enumerate(ii[j]):
                            F[_i][comp] += ff[j][ik]
                    else:
                        F[ii[j], comp] += ff[j]
        elif ind == 5:
            i1, i2, i3 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int)
            p1, p2, p3 = pos[i1], pos[i2], pos[i3]
            k, a0 = par[:,4], par[:,5]
           # a02 = np.full((k.size), np.pi*.4)
            a02 = np.full((k.size), np.pi*.5)
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            ii = [i1, i2, i3]
            _f = - np.heaviside(angle-a02, 1) * k * (angle - a0) * (-1./np.sqrt(1.-(arg)**2))
            for comp in range(2):
                f1 = _f * ( (p2[:,comp]-p3[:,comp])/n_p12/n_p23 + (p2[:,comp]-p1[:,comp])*arg/n_p12**2 )
                f2 = _f * ( (p1[:,comp]-2*p2[:,comp]+p3[:,comp])/n_p12/n_p23 + (p3[:,comp]-p2[:,comp])*arg/n_p23**2 -
                        (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
                f3 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p23 - (p3[:,comp]-p2[:,comp])*arg/n_p23**2 )
                ff = [f1, f2, f3]
                for j in range(3):
                    if len(set(ii[j])) < ii[j].size:
                        for ik, _i in enumerate(ii[j]):
                            F[_i][comp] += ff[j][ik]
                    else:
                        F[ii[j], comp] += ff[j]
        elif ind == 6:
            i1, i2, i3, i4 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int), par[:,4].astype(int)
            p1, p2, p3, p4 = pos[i1], pos[i2], pos[i3], pos[i4]
            k, a0 = par[:,5], par[:,6]
            n_p12, n_p34 = norm(p1-p2, axis=1), norm(p4-p3, axis=1)
            arg = np.multiply(p2-p1, p4-p3).sum(1)/n_p12/n_p34
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            _f = - k * (angle - a0) * (-1./np.sqrt(1.-(arg)**2))
            ii = [i1, i2, i3, i4]
            for comp in range(2):
          #      f1 = _f * ( (p2[:,comp]-p3[:,comp])/n_p12/n_p23 + (p2[:,comp]-p1[:,comp])*arg/n_p12**2 )
          #      f2 = _f * ( (p1[:,comp]-2*p2[:,comp]+p3[:,comp])/n_p12/n_p23 + (p3[:,comp]-p2[:,comp])*arg/n_p23**2 -
          #              (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
          #      f3 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p23 - (p3[:,comp]-p2[:,comp])*arg/n_p23**2 )

                f1 = _f * ( (p3[:,comp]-p4[:,comp])/n_p12/n_p34 + (p2[:,comp]-p1[:,comp])*arg/n_p12**2 )
                f2 = _f * ( ((p1[:,comp]-p2[:,comp])-(p3[:,comp]-p4[:,comp]))/n_p12/n_p34 + (p4[:,comp]-p3[:,comp])*arg/n_p34**2 -
                        (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
            #    f3 = _f * ( (-(p1[:,comp]-p2[:,comp])+(p3[:,comp]-p4[:,comp]))/n_p12/n_p34 + (p4[:,comp]-p3[:,comp])*arg/n_p34**2 -
            #            (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
                f3 = -f1
                f4 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p34 - (p4[:,comp]-p3[:,comp])*arg/n_p34**2 )
                ff = [f1, f2, f3, f4]
                for j in range(4):
                    if len(set(ii[j])) < ii[j].size:
                        for ik, _i in enumerate(ii[j]):
                            F[_i][comp] += ff[j][ik]
                    else:
                        F[ii[j], comp] += ff[j]
        elif ind == 7:
            i1, i2 = __i1, __i2
            p1, p2 = pos[i1], pos[i2]
            p12 = p1-p2
            k = __k_rep_lr
            n_p12 = norm(p12, axis=1)
            f = np.zeros((n, n, 2))
            n12 = np.column_stack((n_p12, n_p12))
            diff = n12 - d
            f[i1, i2] =  k * np.power(n12, -3) * p12
            F += np.sum(f, axis=1)
            F -= np.sum(f, axis=0)
        else:
            print("Don't have potential %d" % ind)
#        print(ind, datetime.datetime.now()-ts)
    return F

