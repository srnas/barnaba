#!/usr/bin/env python

import re, sys, math

import numpy as np
from numpy.linalg import *
from numpy.random import *
from sec_str_constants import *



def energy(pos, param, write_energy=False):
	E = 0. 

	for p in param:
		_type = p[0]
		p1 = pos[p[1]]
		p2 = pos[p[2]]
		n_p12 = norm(p1-p2)
		if _type == 0 or (_type == 1 and n_p12 < d):
			k = p[3]
			d = p[4]
			if k == 0:
				continue
			E_i = 0.5 * k * (n_p12 - d)**2
			E += E_i
		elif _type == 2:
			p3 = pos[p[3]]
			k = p[4]
			if k == 0:
				continue
			a0 = p[5]
			arg = np.dot(p2-p1, p3-p2)/norm(p2-p1)/norm(p3-p2)
			angle = np.arccos(arg)
			E_i = 0.5 * k * (angle - a0)**2
			E += E_i
		else:
			continue
		if write_energy:
			print p, E_i
		
	return E

def force(pos, param, write_force):
	import datetime
	start = datetime.datetime.now()
	_shape = pos.shape
	F = np.zeros(_shape)
	for p in param:
		_type = p[0]
		i1 = p[1]
		i2 = p[2]
		p1 = pos[i1]
		p2 = pos[i2]
		if _type in [0, 1]:
			k = p[3]
			if k == 0:
				continue
			d = p[4]
			n_p12 = norm(p1-p2)
			if _type == 0 or (_type == 1 and n_p12 < d):
				for comp in range(2):
					f = - k * (n_p12 - d) * (p1[comp]-p2[comp]) / n_p12
					F[i1][comp] += f
					F[i2][comp] -= f
					if write_force:
						if _type == 0:
							if k == k_stack1:
								print "Seq_Stacking_F: ",
								print f
							elif k == k_seq:
								print "Sequencial_F: ",
								print f
							elif k == k_bp:
								print "Non-canonical_BP_F: ",
								print f
							else:
								print "WC_BP_F: ",
								print f
								
								
						elif _type == 1:
							print "Repulsive_F(%d_%d): " % (i1, i2),
							print f
						
		elif _type == 2:
			i3 = p[3]
			p3 = pos[i3]
			k = p[4]
			if k == 0:
				continue
			a0 = p[5]
			n_p12 = norm(p2-p1)
			n_p23 = norm(p3-p2)
			arg = np.dot(p2-p1, p3-p2)/n_p12/n_p23
			angle = np.arccos(arg)
			_f = - k * (angle - a0) * (-1./np.sqrt(1-(arg)**2))
			
	#		print "%d %d %d %.0f %.2e" % (i1, i2, i3, angle/np.pi*180,  _f)
			for comp in range(2):
				f1 = _f * ( (p2[comp]-p3[comp])/n_p12/n_p23 + (p2[comp]-p1[comp])*arg/n_p12**2 )
				F[i1][comp] += f1
				f2 = _f * ( (p1[comp]-2*p2[comp]+p3[comp])/n_p12/n_p23 + (p3[comp]-p2[comp])*arg/n_p23**2 -
					(p2[comp]-p1[comp])*arg/n_p12**2 ) 
				F[i2][comp] += f2
				f3 = _f * ( (p2[comp]-p1[comp])/n_p12/n_p23 - (p3[comp]-p2[comp])*arg/n_p23**2 )
				F[i3][comp] += f3
				if write_force:
					print "Angle_F1:", f1, 
					if abs(f1) > 60000: 
						print "%d %d %d %.1f" % (i1, i2, i3, angle/np.pi*180)
					else:
						print 
					print "Angle_F2:", f2, 
					if abs(f2) > 60000: 
						print "%d %d %d %.1f" % (i1, i2, i3, angle/np.pi*180)
					else:
						print 
					print "Angle_F3:", f3, 
					if abs(f3) > 60000:
						print "%d %d %d %.1f" % (i1, i2, i3, angle/np.pi*180)
						
					else:
						print 
			
	end = datetime.datetime.now()
	return F

