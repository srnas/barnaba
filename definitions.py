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



pairings = ['WC','WW','WS','WH','HH','HS','HW','SS','SH','SW',"GU"]
op = ['(','[','{','<']
cl = [')',']','}','>']          

known_abbrev = ["A","C","G","U","N","Y","R","%"]

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
term_5prime = ["OP1","OP2","P","OP3"]

# electronegativity                                                                                                 
chi_O = 1.3
chi_C = 0.4
chi_N = 0.85
chi_P = -0.05

# first index is label, second atoms, third karplus parameters
# you can add more to this list if you wish

                                             # nu from Davies, BD Conformations of nucleosides and nucleotides
                                             # Progress in NMR spectroscopy, 1977
j3 = [["H1H2",   ["H1'","C1'","C2'","1H2'"], [10.2,-0.8,0.0]],\
      ["H2H3",   ["1H2'","C2'","C3'","H3'"], [10.2,-0.8,0.0]],\
      ["H3H4",   ["H3'","C3'","C4'","H4'" ], [10.2,-0.8,0.0]],\
                                             # Generalised Karplus equation from Hasnoot, Altona. Thetraedron, 1980
                                             # P1    P2    --   P3    P4   P5   (P6 = 0 with 3 substituents)
      ["H4H5'",  ["H4'","C4'","C5'","1H5'"], [13.22,-0.99,0.0,0.87,-2.46,19.9] ,[chi_C,0.0,chi_O,chi_O]],\
      ["H4H5''", ["H4'","C4'","C5'","2H5'"], [13.22,-0.99,0.0,0.87,-2.46,19.9] ,[chi_C,chi_O,chi_O,0.0]],\
                                            #  HCOP from Lankhorst, Altona, 1984 
      ["1H5P",   ["1H5'","C5'","O5'","P"],  [15.3,-6.1,1.6]],\
      ["2H5P",   ["2H5'","C5'","O5'","P"],  [15.3,-6.1,1.6]],\
      ["H3P",    ["H3'","C3'","O3'","P"],   [15.3,-6.1,1.6]]]

# Rna only for the moment being
heavy_atoms = ["P","OP1","OP2","O1P","O2P",\
                   "O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'",\
                   "N9","C8","N7","C6","N6","O6","C5",\
                   "C4","N4","O4","N3",\
                   "O2","N2","C2","N1" ]

crazy_vec = [-9999.,-9999.-9999.]
# max bond lenght squared (1.8^2)
maxbond_sq = 3.24
# min bond lenght squared (1.1^2)
minbond_sq = 1.21

pyr = ["C","U","dC","dU","dT"]
pur = ["A","G","dA","dG"]
rna = ["A","U","C","G"]
