import numpy as np

list_bp = ["WC", "GU", "WW", "WH", "HW", "WS", "SW", "SS", "SH", "HS", "HH"]
list_bp_ct = list(np.reshape([[i+"c", i+"t"] for i in list_bp], len(list_bp)*2))
list_wc_pairs = ["WCc", "GUc", "WWc"]
list_stackings = [">>", "<<", "<>", "><"]
list_ann = list_stackings + list_bp_ct

fac = 1
#k_bp = 150 * fac 
k_bp = 150 * fac 
k_wc = 300 * fac
# horizontal interaction in stem
d_short = 32.5

# n, n+2 stacking 
d_stack = 20.
#k_stack = 150 
k_stack = 300 * fac

# n, n+1
d_seq = 20.
k_seq = 2000  
k_seq = 3000  

# diagonal interactions in stem
d_long = np.sqrt(d_short**2+d_seq**2)

#	k_stack = 0

angle = np.pi * .5
k_ang =  10000000

#k_angle_straight = 600000 
k_angle_straight = 1000000 
#k_angle_straight_end = 700000 
k_angle_straight_end = 20000 

k_angle_bulge = 20000 
#k_angle_bulge_rep = 300000 
k_angle_bulge_rep = 100000 
#k_angle_bulge = 200000 


h = .7

#d_rep2 = d_seq * .8
d_rep2 = d_seq 
#d_rep2 = d_seq * .5
#k_rep2 = 5000  
k_rep2 = 10000  
k_rep_lr = 5e4 
k_rep_lr_end = 1e4

k_vertical = 0
k_parallel = 0000000

