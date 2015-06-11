f_factors = [5.,5.,3.]
scale = [1./f_factors[0],1./f_factors[1],1./f_factors[2]]

# mean values and covariance matrix for wc-pair calculation
# extracted from empirical distribution
wc_mean = [2.86,4.67,0.01]
wc_sigma = [[ 0.26, -0.14,0.0],\
            [-0.14,0.13,0.0],\
            [ 0.0, 0.0,  0.33 ]]
det_sigma=0.004686
inv_sigma=[[  9.15492958,   9.85915493,   0.        ],\
               [  9.85915493,  18.30985915,   0.        ],\
               [  0.        ,   0.         ,  3.03030303]]

# treshold values for base pair edges were obtained
# from the angular distribution

theta1 = 0.16
theta2 = 2.0
theta3 = -2.0




interactions = ['..','>>','<<','<>','><','WC','WW','WS','WH','HH','HS','HW','SS','SH','SW','XX']
pairings = ['WC','WW','WS','WH','HH','HS','HW','SS','SH','SW']
op = ['(','[','{','<']
cl = [')',']','}','>']          


#Tolerance to identify zero-eigenvalues modes, 
#corresponding to translational and rotational degrees of freedom.
#Set it to zero or a negative value to print all the eigenvalues/eigenvectors

tol=0.000001 


# list of all RNA atoms for lcs
rna_lcs = ["C2","C4","C6"]
rna_pucker = ["C4'","O4'","C1'","C2'","C3'"]
rna_torsion = ["P","O5'","C5'","C4'","C3'","O3'","O4'","C1'","N9","C4","N1","C2"]
rna_backbone = ["P","O5'","C5'","C4'","C3'","O3'"]
rna_chi_pur = ["O4'","C1'","N9","C4"]
rna_chi_pyr = ["O4'","C1'","N1","C2"]

# Rna only for the moment being
heavy_atoms = ["P","OP1","OP2","O1P","O2P",\
                   "O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                   "N9","C8","N7","C6","N6","O6","C5",\
                   "C4","N4","O4","N3",\
                   "O2","N2","C2","N1" ]

# max bond lenght squared in sugar (1.75^2)
maxbond_pucker_sq = 3.0625
# min bond lenght squared in sugar (1.1^2)
minbond_pucker_sq = 1.21
