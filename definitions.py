


purines = ["A","G"]
pyrimidines = ["C","U","T"]

residue_dict = {'A': 'A', 'rA':'A','RA':'A','RA5':'A','RA3':'A','A3':'A','A5':'A',\
                'C': 'C', 'rC':'C','RC':'C','RC5':'C','RC3':'C','C3':'C','C5':'C',\
                'G': 'G', 'rG':'G','RG':'G','RG5':'G','RG3':'G','G3':'G','G5':'G',\
                'U': 'U', 'rU':'U','RU':'U','RU5':'U','RU3':'U','U3':'U','U5':'U',\
                'T': 'dT', 'dT':'dT','DT':'dT','DT5':'dT','DT3':'dT',\
                'dA':'dA','DA':'dA','DA5':'dA','DA3':'dA',\
                'dC':'dC','DC':'dC','DC5':'dC','DC3':'dC',\
                'dG':'dG','DG':'dG','DG5':'dG','DG3':'dG'}

modified_dict = {'1MA':'A','5AA':'A','P5P':'A','2MA':'A',\
                'OMC':'C','5MC':'C','CBV':'C',\
                '2MG':'G','YG':'G','7MG':'G','OMG':'G','1MG':'G','M2G':'G',\
                'H2U':'U','PSU':'U','OMU':'U','UR3':'U','5MU':'U','5BU':'U','4SU':'U'}


others = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HSD","HSE","HSP","ILE","HIS","LEU","LYS","LSN","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HOH","MG","K","NA","CL","CD","CA"]


# geometric definitions
f_factors = [0.5,0.5,0.3]
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
complementary = {'A':'U','C':'G','U':'A','G':'C'}
bb_atoms = ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'","P"]
