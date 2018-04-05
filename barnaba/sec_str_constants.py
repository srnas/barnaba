import numpy as np

list_bp = ["WC", "GU", "WW", "WH", "HW", "WS", "SW", "SS", "SH", "HS", "HH"]
list_bp_ct = list(np.reshape([[i+"c", i+"t"] for i in list_bp], len(list_bp)*2))
list_wc_pairs = ["WCc", "GUc", "WWc"]
list_stackings = [">>", "<<", "<>", "><"]
list_ann = list_stackings + list_bp_ct

threshold = 0.1

factor = 100
k_bp = 100 * factor
k_wc = 300 *factor
# horizontal interaction in stem
d_short = 32.5

# n, n+2 stacking 
d_stack = 20.
k_stack = 150 * factor

# n, n+1
d_seq = 20.
k_seq = 2000 * factor 

# diagonal interactions in stem
d_long = np.sqrt(d_short**2+d_seq**2)

#	k_stack = 0

angle = np.pi * .5
k_ang = 200000 * factor

k_angle_180 = 300000 * factor
#k_vertical = 100000 * factor
k_vertical = 200000 * factor



#d_rep2 = d_seq * .8
d_rep2 = d_seq * .8
k_rep2 = 5000 * factor 

k_rep1 = k_rep2*.5
k_rep3 = k_rep1*.01
k_angle_180 = 0

d_pull = 10*d_seq
k_pull = 0 * factor

F_max = .1 * factor

