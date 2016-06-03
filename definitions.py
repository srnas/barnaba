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

known_abbrev = ["A","C","G","U","N","Y","R","%"]

#Tolerance to identify zero-eigenvalues modes, 
#corresponding to translational and rotational degrees of freedom.
#Set it to zero or a negative value to print all the eigenvalues/eigenvectors

tol=0.000001 


# list of all RNA atoms for lcs
rna_lcs = ["C2","C4","C6"]
rna_pucker = ["C4'","O4'","C1'","C2'","C3'"]
rna_torsion_pur = ["P","O5'","C5'","C4'","C3'","O3'","O4'","C1'","C2'","N9","C4"]
rna_torsion_pyr = ["P","O5'","C5'","C4'","C3'","O3'","O4'","C1'","C2'","N1","C2"]

#rna_backbone = ["P","O5'","C5'","C4'","C3'","O3'"]


# electronegativity                                                                                                 
chi_O = 1.3
chi_C = 0.4
chi_N = 0.85
chi_P = -0.05

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

# Rna only for the moment being
heavy_atoms = ["P","OP1","OP2","O1P","O2P",\
                   "\"O5'\"","\"C5'\"","\"C4'\"","\"O4'\"","\"C3'\"","\"O3'\"","\"C2'\"","\"O2'\"","\"C1'\"",\
                   "N9","C8","N7","C6","N6","O6","C5",\
                   "C4","N4","O4","N3",\
                   "O2","N2","C2","N1" ]

align_atoms_pur = ["\"O5'\"","\"C5'\"","\"C4'\"","\"O4'\"","\"C3'\"","\"O3'\"","\"C2'\"","\"O2'\"","\"C1'\"",\
                   "C4","C5","C6","N1","C2","N3"]
align_atoms_pyr = ["\"O5'\"","\"C5'\"","\"C4'\"","\"O4'\"","\"C3'\"","\"O3'\"","\"C2'\"","\"O2'\"","\"C1'\"",\
                   "N1","C6","C5","C4","N3","C2"]

#crazy_vec = [-9999.,-9999.-9999.]
# max bond lenght squared (1.8^2)
#maxbond_sq = 3.24
# min bond lenght squared (1.1^2)
#minbond_sq = 1.21

pyr = ["C","U","rU","RU","RU5","RU3","U3","U5","rC","RC","RC5","RC3","C3","C5"]
pur = ["A","G","rG","RG","RG5","RG3","G3","G5","rA","RA","RA5","RA3","A3","A5"]

# if you add here, add also to pypu_dict and residue_dict
rna = ["U","rU","RU","RU5","RU3","U3","U5",\
       "C","rC","RC","RC5","RC3","C3","C5",\
       "G","rG","RG","RG5","RG3","G3","G5",\
       "A","rA","RA","RA5","RA3","A3","A5",\
       "\"1MA\"","OMC","\"5MC\"","\"2MG\"","M2G","OMG","YG",\
       "PSU","H2U","\"5MU\""]


complementary = {'A':'U','C':'G','U':'A','G':'C'}

pypu_dict = {'U': 'Y', 'rU':'Y','RU':'Y','RU5':'Y','RU3':'Y','U3':'Y','U5':'Y',\
                'H2U':'Y','PSU':'Y','OMU':'Y','UR3':'Y','5MU':'Y',\
                'A': 'R', 'rA':'R','RA':'R','RA5':'R','RA3':'R','A3':'R','A5':'R',\
                '1MA':'R',\
                'C': 'Y', 'rC':'Y','RC':'Y','RC5':'Y','RC3':'Y','C3':'Y','C5':'Y',\
                'OMC':'Y','5MC':'Y',\
                'G': 'R', 'rG':'R','RG':'R','RG5':'R','RG3':'R','G3':'R','G5':'R',\
                'M2G':"R",'2MG':'R','YG':'R','7MG':'R','OMG':'R','1MG':'R'}
            
residue_dict = {'U': 'U', 'rU':'U','RU':'U','RU5':'U','RU3':'U','U3':'U','U5':'U',\
                'H2U':'U','PSU':'U','OMU':'U','UR3':'U','5MU':'U',\
                'A': 'A', 'rA':'A','RA':'A','RA5':'A','RA3':'A','A3':'A','A5':'A',\
                '1MA':'A',\
                'C': 'C', 'rC':'C','RC':'C','RC5':'C','RC3':'C','C3':'C','C5':'C',\
                'OMC':'C','5MC':'C',\
                'G': 'G', 'rG':'G','RG':'G','RG5':'G','RG3':'G','G3':'G','G5':'G',\
                '2MG':'G','YG':'G','7MG':'G','OMG':'G','1MG':'G','M2G':'G',\
                'T': 'dT', 'dT':'dT','DT':'dT','DT5':'dT','DT3':'dT',\
                'dA':'dA','DA':'dA','DA5':'dA','DA3':'dA',\
                'dC':'dC','DC':'dC','DC5':'dC','DC3':'dC',\
                'dG':'dG','DG':'dG','DG5':'dG','DG3':'dG',\
                "ALA":"ALA","ARG":"ARG","ASN":"ASN","ASP":"ASP",\
                "CYS":"CYS","GLN":"GLN","GLU":"GLU","GLY":"GLY",\
                "HSD":"HSD","HSE":"HSD","HSP":"HSD","ILE":"ILE","HIS":"HIS",\
                "LEU":"LEU","LYS":"LYS","LSN":"LYS","MET":"MET",\
                "PHE":"PHE","PRO":"PRO","SER":"SER","THR":"THR",\
                "TRP":"TRP","TYR":"TYR","VAL":"VAL"}



annotation = ["..","XX","WWt","WHt","WSt","HHt","HWt","HSt","SHt","SWt","SSt","WCt",\
              "WWc","WHc","WSc","HHc","HWc","HSc","SHc","SWc","SSc","WCc",\
              ">>","<<","><","<>"]

aa = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HSD","HSE","HSP","ILE","HIS","LEU","LYS","LSN","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]


#"HOH":"HOH","MG":"MG","K":"K","NA":"NA","CL":"CL","CD":"CD","CA":"CA",\
#    "BRU":"BRU","RHD":"RHD","GDP":"GDP","A23":"A23","B12":"B12","NME":"NME",\
#"PGP":"PGP","LU":"LU","SO4":"SO4","M2G":"M2G","0C":"0C","0G":"0G","0U":"0U",\
#                "3CO":"3CO","42B":"42B","5BU":"5BU","5IC":"5IC","A2M":"A2M","A5M":"A5M","AM2":"AM2",\
#                "BA":"BA","BDG":"BDG","BGM":"BGM","BR":"BR","BTN":"BTN","CBR":"CBR","CBV":"CBV",\
#                "CNC":"CNC","CO":"CO","CSL":"CSL","DAI":"DAI","GET":"GET","0C":"0C","0G":"0G","0U":\
#                "0U","3CO":"3CO","42B":"42B","5BU":"5BU","5IC":"5IC","A2M":"A2M","A5M":"A5M",\
#                "AM2":"AM2","BA":"BA","BDG":"BDG","BGM":"BGM","BR":"BR","BTN":"BTN","CBR":"CBR",\
#                "CBV":"CBV","CNC":"CNC","CO":"CO","CSL":"CSL","DAI":"DAI","GET":"GET","GTP":"GTP",\
#                "I":"I","MN":"MN","MPD":"MPD","MTU":"MTU","NCO":"NCO","NH4":"NH4","NMY":"NMY",\
#                "P24":"P24","PA2":"PA2","PAR":"PAR","PB":"PB","PO2":"PO2","PPU":"PPU","PT4":"PT4",\
#                "PYY":"PYY","RIA":"RIA","ROS":"ROS","SPM":"SPM","SR":"SR","SRY":"SRY","T6A":"T6A",\
#                "TL":"TL","TOY":"TOY","UMS":"UMS","YYG":"YYG","ZN":"ZN","2AD":"2AD","3AY":"3AY",\
#                "6AP":"6AP","AB6":"AB6","AB9":"AB9","ACT":"ACT","AKN":"AKN","CNY":"CNY","HPA":"HPA",\
#                "JS4":"JS4","JS5":"JS5","KAN":"KAN","LIV":"LIV","LLL":"LLL","N6G":"N6G","NF2":"NF2",\
#                "OS":"OS","P1P":"P1P","RIO":"RIO","S4C":"S4C","TPP":"TPP","XXX":"XXX","MES":"MES",\
#                "CCC":"CCC","IRI":"IRI","SAM":"SAM","U33":"U33","1PE":"1PE","1SC":"1SC","2AU":"2AU","2BP":"2BP","3AD":"3AD","3AW":"3AW","3DA":"3DA","5AZ":"5AZ","5CF":"5CF","5GP":"5GP","6GO":"6GO","6GU":"6GU","A2F":"A2F","A6A":"A6A","A6C":"A6C","A6G":"A6G","A6U":"A6U","AF2":"AF2","AT7":"AT7","AU3":"AU3","BFT":"BFT","C2E":"C2E","CAC":"CAC","CFZ":"CFZ","CS":"CS","D2X":"D2X","DGP":"DGP","DX4":"DX4","EEM":"EEM","F":"F","FFO":"FFO","FMN":"FMN","FOZ":"FOZ","G46":"G46","G6P":"G6P","GF2":"GF2","GLP":"GLP","GMP":"GMP","GNG":"GNG","GOL":"GOL","GRB":"GRB","HRG":"HRG","IEL":"IEL","IOD":"IOD","IPA":"IPA","JS6":"JS6","LCA":"LCA","LCC":"LCC","LCG":"LCG","LHA":"LHA","N30":"N30","N33":"N33","OHX":"OHX","OLZ":"OLZ","PO4":"PO4","PQ0":"PQ0","PRF":"PRF","PYI":"PYI","R14":"R14","RS3":"RS3","RUS":"RUS","S9L":"S9L","SAH":"SAH","SE4":"SE4","SFG":"SFG","SIN":"SIN","SLZ":"SLZ","SS0":"SS0","TB":"TB","THF":"THF","TLN":"TLN","TPS":"TPS","UFT":"UFT","US5":"US5","UZR":"UZR","XAN":"XAN","XUG":"XUG","6HS":"6HS","A44":"A44","ACA":"ACA","C43":"C43","EPE":"EPE","G48":"G48","I2A":"I2A","IR3":"IR3","IU":"IU","P5P":"P5P","SIS":"SIS","U36":"U36","MGT":"MGT","4SU":"4SU","U34":"U34"};
#
