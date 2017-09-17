import re
import itertools as its

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

known_abbrev = ["A","C","G","U","N","Y","R","%"]
term_5prime = ["OP1","OP2","P","OP3"]

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
bb_atoms = ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'","P","OP1","OP2"]
nt_atoms = {"A":bb_atoms + ["N1","C2","N3","C4","C5","C6","N6","N7","C8","N9"],\
            "G":bb_atoms + ["N1","C2","N2","N3","C4","C5","C6","O6","N7","C8","N9"],\
            "C":bb_atoms + ["N1","C2","O2","N3","C4","N4","C5","C6"],\
            "U":bb_atoms + ["N1","C2","O2","N3","C4","O4","C5","C6"]}
rna = ["A","C","G","U","T"]

# Scalar couplitg
# first index is label, second atoms, third karplus parameters
# you can add more to this list if you wish
td_pi = 2.094395
# nu from Davies, BD Conformations of nucleosides and nucleotides
# Progress in NMR spectroscopy, 1977
j3 = [["H1H2",   ["H1'","C1'","C2'","1H2'"], [10.2,-0.8,0.0,0.0,0.0]],\
      ["H2H3",   ["1H2'","C2'","C3'","H3'"], [10.2,-0.8,0.0,0.0,0.0]],\
      ["H3H4",   ["H3'","C3'","C4'","H4'" ], [10.2,-0.8,0.0,0.0,0.0]],\
      # Generalised Karplus equation from Hasnoot, Altona. Thetraedron, 1980
      # P1    P2    --   P3    P4   P5   (P6 = 0 with 3 substituents)
      #["H4H5'",  ["H4'","C4'","C5'","1H5'"], [13.22,-0.99,0.0,0.87,-2.46,19.9] ,[chi_C,0.0,chi_O,chi_O]],\
      #["H4H5''", ["H4'","C4'","C5'","2H5'"], [13.22,-0.99,0.0,0.87,-2.46,19.9] ,[chi_C,chi_O,chi_O,0.0]],\

      # simiplification of altona hasnoot - in the form  A*cos*cos + B*cos + C + D*sin*cos
      ["H4H5'",  ["C3'","C4'","C5'","O5'"], [8.313139, -0.99, 1.373430,0.269906,-td_pi]],\
      ["H4H5''", ["C3'","C4'","C5'","O5'"], [8.313139, -0.99, 1.373430,-4.752290,0.0]],\
        #  HCOP from Lankhorst, Altona, 1984 
      #["1H5P",   ["1H5'","C5'","O5'","P"],  [15.3,-6.1,1.6]],\
      #["2H5P",   ["2H5'","C5'","O5'","P"],  [15.3,-6.1,1.6]],\
      #  HCOP from Lee, Sarma 1976 
      ["1H5P",   ["C4'","C5'","O5'","P"],  [18.1,-4.8,1.5,0.0,-td_pi]],\
      ["2H5P",   ["C4'","C5'","O5'","P"],  [18.1,-4.8,1.5,0.0,td_pi]],\
      # Marino and Scwhalbe
      ["C4Pb",    ["C4'","C5'","O5'","P"],  [6.9,-3.4,0.7,0.0,0.0]],\
      ["C4Pe",    ["C4'","C3'","O3'","P"],  [6.9,-3.4,0.7,0.0,0.0]],\
      #  HCOP from Lankhorst, Altona, 1984 
      ["H3P",    ["C4'","C3'","O3'","P"],   [15.3,-6.1,1.6,0.0,td_pi]]]


def get_pattern(query):
    # build pattern for regular expression
    pattern = "^"
    for res in query:
        assert res in known_abbrev, "# Fatal error: character %s not known. Use AUCG/NYR" % (res)
        if(res in rna):
            pattern += res
        else:
            if(res == "N"):
                pattern += "[AUCGT]"
            if(res == "Y"):
                pattern += "[UCT]"
            if(res == "R"):
                pattern += "[AG]"
    pattern += "$"
    return pattern

def get_idx(sequence,query,bulges=0):

    ll = len(query)
    seq_str  = "".join(sequence)
    pattern = get_pattern(query)

    # generate first all possible indeces
    indeces = []
    for b in range(bulges+1):

        # create position of insertions
        comb= its.combinations(range(1,ll+b-1),b)
        for it1 in comb:
            idx1 = range(ll+b)
            # remove them from list
            for it2 in it1:
                idx1.remove(it2)
            for start in range(len(sequence)-ll-b+1):
                # get sequence
                idx2 =  [idx1[i] + start for i in range(ll)]
                substr = "".join([sequence[i] for i in idx2])
                if(re.match(pattern,substr) != None):
                    indeces.append(idx2)
    return indeces
