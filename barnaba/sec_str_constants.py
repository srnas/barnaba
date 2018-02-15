import numpy as np

list_bp = ["WC", "GU", "WW", "WH", "HW", "WS", "SW", "SS", "SH", "HS", "HH"]
list_bp_ct = list(np.reshape([[i+"c", i+"t"] for i in list_bp], len(list_bp)*2))
list_wc_pairs = ["WCc", "GUc", "WWc"]
list_stackings = [">>", "<<", "<>", "><"]
list_ann = list_stackings + list_bp_ct

threshold = 0.1

factor = 100
d_bp = 32.5
#k_bp = 300
#k_bp = 500
k_bp = 200 * factor
#k_bp = 2000
#	k_bp = 0

# sequential stacking
d_stack1 = 20
k_stack1 = 0

d_seq = 20
#	k_seq = 20000 
k_seq = 200 * factor 

# non-sequential stacking
#d_stack2 = d_bp
d_stack2 = np.sqrt(d_bp**2+d_seq**2)
# lonely basepair
d_bp2 = d_stack2
k_stack2 = k_bp 

#	k_stack = 0

angle = np.pi * .5
k_ang = 5000 * factor



d_rep1 = np.sqrt(d_bp**2+d_seq**2)
#d_rep = d_seq
#	k_rep = 10000
k_rep1 = k_bp
#k_rep1 = 100 * factor


d_rep2 = d_seq * .8
#k_rep2 = 500 * factor
k_rep2 = 200 * factor 
	
F_max = .1 * factor

