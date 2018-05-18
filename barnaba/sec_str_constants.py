import numpy as np

list_bp = ["WC", "GU", "WW", "WH", "HW", "WS", "SW", "SS", "SH", "HS", "HH"]
list_bp_ct = list(np.reshape([[i+"c", i+"t"] for i in list_bp], len(list_bp)*2))
list_wc_pairs = ["WCc", "GUc", "WWc"]
list_stackings = [">>", "<<", "<>", "><"]
list_ann = list_stackings + list_bp_ct

k_bp = 150 
k_wc = 300
# horizontal interaction in stem
d_short = 32.5

# n, n+2 stacking 
d_stack = 20.
k_stack = 300 

# n, n+1
d_seq = 20.
k_seq = 3000  

# diagonal interactions in stem
d_long = np.sqrt(d_short**2+d_seq**2)

#	k_stack = 0

angle = np.pi * .5
k_ang =  1000
k_ang_end =  500000

k_angle_straight =  1500000 
k_angle_straight_end =   10000 

h = 1.

d_rep2 = d_seq 
k_rep2 = 10000  
k_rep_lr =  5e4 
k_rep_lr_end = 0.# 1e4


