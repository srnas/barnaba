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
            e_i = .5 * k * (n_p12 - d)**2
            E_i += sum(e_i)
        elif ind == 1:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k, d = par[:,3], par[:,4]
            n_p12 = norm(p1-p2, axis=1)
            e_i = np.heaviside(d-n_p12, 1) * .5 * k * (n_p12 - d)**2
            E_i += sum(e_i)
        elif ind == 2:
            i1, i2, i3 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int)
            p1, p2, p3 = pos[i1], pos[i2], pos[i3]
            k, a0 = par[:,4], par[:,5]
            n_p12, n_p23 = norm(p1-p2, axis=1), norm(p3-p2, axis=1)
            arg = np.multiply(p2-p1, p3-p2).sum(1)/n_p12/n_p23
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            e_i = 0.5 * k * (angle - a0)**2
            E_i += sum(e_i)
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
            e_i = 0.5 * k * (angle - a0)**2
            E_i += sum(e_i)
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
            e_i = np.heaviside(a02-angle, 1) * 0.5 * k * (angle - a0)**2
            E_i += sum(e_i)
        elif ind == 6:
            i1, i2, i3, i4 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int), par[:,4].astype(int)
            p1, p2, p3, p4 = pos[i1], pos[i2], pos[i3], pos[i4]
            k, a0 = par[:,5], par[:,6]
            n_p12, n_p34 = norm(p1-p2, axis=1), norm(p4-p3, axis=1)
            arg = np.multiply(p2-p1, p4-p3).sum(1)/n_p12/n_p34
            arg[np.where(arg >= 1.)] = 1.-1e-6 
            arg[np.where(arg <= -1.)] = -1.+1e-6 
            angle = np.arccos(arg)
            e_i = 0.5 * k * (angle - a0)**2
            E_i += sum(e_i)

    return E_i

def force(pos, param, write_force):
    _shape = pos.shape
    F = np.zeros(_shape)
    for ind, par in param.items():
        if ind == 0:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k, d = par[:,3], par[:,4]
            n_p12 = norm(p1-p2, axis=1)
            for comp in range(2):
                f = - k * (n_p12 - d) * (p1[:,comp]-p2[:,comp]) / n_p12
                if len(set(i1)) < i1.size:
                    for ik, _i1 in enumerate(i1):
                        F[_i1][comp] += f[ik]
                else:
                    F[i1, comp] += f
                if len(set(i2)) < i2.size:
                    for ik, _i2 in enumerate(i2):
                        F[_i2][comp] -= f[ik]
                else:
                    F[i2, comp] -= f

        elif ind == 1:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
            p1, p2 = pos[i1], pos[i2]
            k, d = par[:,3], par[:,4]
            n_p12 = norm(p1-p2, axis=1)
            for comp in range(2):
                f = - np.heaviside(d-n_p12, 1) * k * (n_p12 - d) * (p1[:,comp]-p2[:,comp]) / n_p12
                if len(set(i1)) < i1.size:
                    for ik, _i1 in enumerate(i1):
                        F[_i1][comp] += f[ik]
                else:
                    F[i1, comp] += f
                if len(set(i2)) < i2.size:
                    for ik, _i2 in enumerate(i2):
                        F[_i2][comp] -= f[ik]
                else:
                    F[i2, comp] += f
        elif ind == 2:
            i1, i2, i3 = par[:,1].astype(int), par[:,2].astype(int), par[:,3].astype(int)
            p1, p2, p3 = pos[i1], pos[i2], pos[i3]
            k, a0 = par[:,4], par[:,5]
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
                f3 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p23 - (p3[:,comp]-p2[:,comp])*arg/n_p23**2 )
                ii = [i1, i2, i3]
                ff = [f1, f2, f3]
                for j in range(3):
                    if len(set(ii[j])) < ii[j].size:
                        for ik, _i in enumerate(ii[j]):
                            F[_i][comp] += ff[j][ik]
                    else:
                        F[ii[j], comp] += ff[j]
        elif ind == 4:
            i1, i2 = par[:,1].astype(int), par[:,2].astype(int)
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
                ii = [i1, i2]
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
            _f = - np.heaviside(angle-a02, 1) * k * (angle - a0) * (-1./np.sqrt(1.-(arg)**2))
            for comp in range(2):
                f1 = _f * ( (p2[:,comp]-p3[:,comp])/n_p12/n_p23 + (p2[:,comp]-p1[:,comp])*arg/n_p12**2 )
                f2 = _f * ( (p1[:,comp]-2*p2[:,comp]+p3[:,comp])/n_p12/n_p23 + (p3[:,comp]-p2[:,comp])*arg/n_p23**2 -
                        (p2[:,comp]-p1[:,comp])*arg/n_p12**2 ) 
                f3 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p23 - (p3[:,comp]-p2[:,comp])*arg/n_p23**2 )
                ii = [i1, i2, i3]
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
                f3 = 0
                f4 = _f * ( (p2[:,comp]-p1[:,comp])/n_p12/n_p34 - (p4[:,comp]-p3[:,comp])*arg/n_p34**2 )
                ii = [i1, i2, i3, i4]
                ff = [f1, f2, f3, f4]
                for j in range(4):
                    if len(set(ii[j])) < ii[j].size:
                        for ik, _i in enumerate(ii[j]):
                            F[_i][comp] += ff[j][ik]
                    else:
                        F[ii[j], comp] += ff[j]
        else:
            print("Don't have potential %d" % ind)
    return F

