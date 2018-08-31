#   This is baRNAba, a tool for analysis of nucleic acid 3d structure
#   Copyright (C) 2017 Sandro Bottaro (sandro.bottaro@bio.ku.dk)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import, division, print_function
import re
import itertools as its
import collections
import numpy as np

tol=1.0e-06 
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

known_abbrev = ["A","C","G","U","N","Y","R","S","W","K","M","B","D","H","V","%"]
term_5prime = ["OP1","OP2","P","OP3"]

# geometric definitions
f_factors = [0.5,0.5,0.3]
scale = [1./f_factors[0],1./f_factors[1],1./f_factors[2]]


# donors

donors = {"A":["N6","C2","C8","O2'"],\
          "C":["N4","C5","C6","O2'"],\
          "G":["N1","N2","C8","O2'"],\
          "U":["N3","C5","C6","O2'"],\
          "dA":["N6","C2","C8"],\
          "dC":["N4","C5","C6"],\
          "dG":["N1","N2","C8"],\
          "dT":["N3","C5","C6"]}

acceptors = {"A":["N1","N3","N7","O2'"],\
             "C":["N3","O2","O2'"],\
             "G":["O6","N3","N7","O2'"],\
             "U":["O2","O4","O2'"],\
             "dA":["N1","N3","N7"],\
             "dC":["N3","O2"],\
             "dG":["O6","N3"],\
             "dT":["O2","O4"]}

glyco = {"A":["C1'","N9"],
         "C":["C1'","N1"],
         "G":["C1'","N9"],
         "U":["C1'","N1"],
         "dA":["C1'","N9"],
         "dC":["C1'","N1"],
         "dG":["C1'","N9"],
         "dT":["C1'","N1"]}

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

bb_angles = ["alpha","beta","gamma","delta","eps","zeta","chi"]
sugar_angles = ["nu1","nu2","nu3","nu4","nu5"]
pucker = ["phi","amp"]
# this mapping serves to get the correct index 
couplings_idx = collections.OrderedDict([("H1H2",0),("H2H3",1),("H3H4",2),\
                             ("1H5P",3),("2H5P",3),("C4Pb",3),\
                             ("1H5H4",4),("2H5H4",4),\
                             ("H3P",5),("C4Pe",5),\
                             ("H1C2/4",6),("H1C6/8",6)])
couplings_karplus = {\
                     
                     ####### SUGAR  ###########
                     # Davies, BD Conformations of nucleosides and nucleotides  Progress in NMR spectroscopy, 1977
                     #"H1H2": [10.2,-0.8,0.0,0.0,0.0],\
                     #"H2H3": [10.2,-0.8,0.0,0.0,0.0],\
                     #"H3H4": [10.2,-0.8,0.0,0.0,0.0],\
                     # Generalised Karplus equation from Hasnoot, Altona. Thetraedron, 1980
                     #"H1H2": [6.96462,-0.91,1.02629,1.27009,0.0],\  
                     #"H2H3": [8.28926,-0.91,0.66772,0.00193297,0.0],\
                     #"H3H4": [7.96446,-0.91,0.77241,-0.262475,0.0],\
                     ## Condon, Stacking in rna: Nmr of four tetramers benchmark molecular dynamics. Journal of Chemical Theory and Computation, 2015.
                     "H1H2": [9.67,-2.03,0.0,0.0,0.0],\
                     "H2H3": [9.67,-2.03,0.0,0.0,0.0],\
                     "H3H4": [9.67,-2.03,0.0,0.0,0.0],\
                     ######### BETA #####
                     #  HCOP from Lee, Sarma 1976 
                     #"1H5P":[18.1,-4.8,1.5,0.0,-td_pi],\
                     #"2H5P":[18.1,-4.8,1.5,0.0,td_pi],\
                     # HCOP from Lankhorst, Journal of Biomolecular Structure and Dynamics, 1(6):1387 1405, 1984.
                     "1H5P":[15.3,-6.1,1.6,0.0,-td_pi],\
                     "2H5P":[15.3,-6.1,1.6,0.0,td_pi],\
                     # HCOP from Mooren Nucleic acids research, 22(13):2658 2666, 1994
                     #"1H5P":[15.3,-6.2,1.5,0.0,-td_pi],\
                     #"2H5P":[15.3,-6.2,1.5,0.0,td_pi],\
                     #HCOP from Marino, Accounts of chemical research, 32(7):614 623, 1999.
                     "C4Pb":[6.9,-3.4,0.7,0.0,0.0],\
                     #HCOP from wijmenga, Progress in nuclear magnetic resonance spectroscopy, 32(4):287 387, 1998
                     #"C4Pb":[8.0,-3.4,0.5,0.0,0.0], \
                     ##### GAMMA  #########
                     # from Davies, BD Conformations of nucleosides and nucleotides  Progress in NMR spectroscopy, 1977
                     "1H5H4":[9.7,-1.8,0.0,0.0,-td_pi],\
                     "2H5H4":[9.7,-1.8,0.0,0.0,0.0],\
                     # from Altona, Hasnoot 
                     #"1H5H4":[8.313139, -0.99, 1.373430,0.269906,-td_pi],\
                     #"2H5H4":[8.313139, -0.99, 1.373430,-4.752290,0.0],\
                     ##### EPSILON ########
                     # From Lankhorst, Altona, 1984
                     "H3P":[15.3,-6.1,1.6,0.0,td_pi],\
                     # Marino and Scwhalbe  Accounts of chemical research, 32(7):614 623, 1999
                     "C4Pe":[6.9,-3.4,0.7,0.0,0.0],\
                     ######  CHI #########
                     # Ippel Magnetic resonance in chemistry, 34(13):S156 S176, 1996.
                     "H1C2/4":[4.7,2.3,0.1,0.0,-0.5*td_pi],\
                     "H1C6/8":[4.5,-0.6,0.1,0.0,-0.5*td_pi],
                     
}

             
# Nomenclature Committee of the International Union of Biochemistry (NC-IUB).
# Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences.
# Recommendations 1984. Biochem. J. 143001985, 229, 281-286.
def get_pattern(query):
    # build pattern for regular expression
    pattern = "^"
    for res in query:
        assert res in known_abbrev, "# Fatal error: character %s not known. Use AUCG/NYRSWKMBDHV" % (res)
        if(res in rna):
            pattern += res
        else:
            if(res == "N"): # aNy
                pattern += "[AUCGT]"
            if(res == "Y"): # pYrimidine
                pattern += "[UCT]"
            if(res == "R"): # puRine
                pattern += "[AG]"
            if(res == "S"): # Strong
                pattern += "[GC]"
            if(res == "W"): # Weak
                pattern += "[AUT]"
            if(res == "K"): # Keto
                pattern += "[GUT]"
            if(res == "M"): # aMino
                pattern += "[AC]"
            if(res == "B"): # not adenine
                pattern += "[UCGT]"
            if(res == "D"): # not cytosine
                pattern += "[AUGT]"
            if(res == "H"): # not guanine
                pattern += "[AUCT]"
            if(res == "V"): # not uracil
                pattern += "[ACG]"

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
            idx1 = list(range(ll+b))
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


def unique_rows(ar):
    """Remove repeated rows from a 2D array.

    In particular, if given an array of coordinates of shape
    (Npoints, Ndim), it will remove repeated points.

    Parameters
    ----------
    ar : 2-D ndarray
        The input array.

    Returns
    -------
    ar_out : 2-D ndarray
        A copy of the input array with repeated rows removed.

    Raises
    ------
    ValueError : if `ar` is not two-dimensional.

    Notes
    -----
    The function will generate a copy of `ar` if it is not
    C-contiguous, which will negatively affect performance for large
    input arrays.

    Examples
    --------
    >>> ar = np.array([[1, 0, 1],
    ...                [0, 1, 0],
    ...                [1, 0, 1]], np.uint8)
    >>> unique_rows(ar)
    array([[0, 1, 0],
           [1, 0, 1]], dtype=uint8)
    """
    if ar.ndim != 2:
        raise ValueError("unique_rows() only makes sense for 2D arrays, "
                         "got %dd" % ar.ndim)
    # the view in the next line only works if the array is C-contiguous
    ar = np.ascontiguousarray(ar)
    # np.unique() finds identical items in a raveled array. To make it
    # see each row as a single item, we create a view of each row as a
    # byte string of length itemsize times number of columns in `ar`
    ar_row_view = ar.view('|S%d' % (ar.itemsize * ar.shape[1]))
    _, unique_row_indices = np.unique(ar_row_view, return_index=True)
    ar_out = ar[unique_row_indices]
    return ar_out
