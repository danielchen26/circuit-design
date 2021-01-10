import numpy as np
import pandas as pd
import scipy as sp
from pyeda.inter import *
from graphviz import Source
import itertools


# Table S1. Parameters for gate response functions
para_s1 = pd.DataFrame()
para_s1['repressor'] = ['AmeR','AmtR','BetI','BM3R1','BM3R1','BM3R1','HlyIIR','IcaRA','LitR','PhIF','PhIF','PhIF','QacR','QacR','SrpR','SrpR','SrpR','SrpR']
para_s1['RBS']       = ['F1','A1','E1','B1','B2','B3','H1','I1','L1','P1','P2','P3','Q1','Q2','S1','S2','S3','S4']
para_s1['Y_min']     = [0.2,0.06,0.07,0.004,0.005,0.01,0.07,0.08,0.07,0.01,0.02,0.02,0.01,0.03,0.003,0.003,0.004,0.007]
para_s1['Y_max']     = [3.8,3.8,3.8,0.5,0.5,0.8,2.5,2.2,4.3,3.9,4.1,6.8,2.4,2.8,1.3,2.1,2.1,2.1]
para_s1['K']         = [0.09,0.07,0.41,0.04,0.15,0.26,0.19,0.1,.05,0.03,0.13,0.23,0.05,0.21,0.01,0.04,0.06,0.1]
para_s1['n']         = [1.4,1.6,2.4,3.4,2.9,3.4,2.6,1.4,1.7,4.0,3.9,4.2,2.7,2.4,2.9,2.6,2.8,2.8]
para_s1.to_csv("./para_s1.csv", index=False)



# Table S4: Insulated gate response function parameters
para_s4 = pd.DataFrame()
para_s4['repressor'] =['AmeR','AmeR','BetI','BM3R1','BM3R1','BM3R1','HlyIIR','IcaRA','LitR','LmrA','PhIF','PhIF','PhIF','PsrA','QacR','QacR','SrpR','SrpR','SrpR','SrpR']
para_s4['RBS'] = ['F1','A1','E1','B1','B2','B3','H1','I1','l1','N1','P1','P2','P3','R1','Q1','Q2','S1','S2','S3','S4']
para_s4['Y_min'] = [0.2,  0.06,  0.07,  0.004,  0.005,  0.01, 0.07, 0.08,  0.07,  0.2,  0.01,  0.02, 0.02,  0.2,  0.01,  0.03, 0.003, 0.003,   0.004, 0.007]
para_s4['Y_max'] = [3.8,  3.8,   3.8,   0.5,    0.5,    0.8,  2.5,  2.2,   4.3,   2.2,  3.9,   4.1,  6.8,   5.9,  2.4,   2.8,  1.3,   2.1,     2.1,   2.1]
para_s4['K'] =     [0.09, 0.07,  0.41,  0.04,   0.15,   0.26, 0.19, 0.10,  0.05,  0.18, 0.03,  0.13, 0.23,  0.19, 0.05,  0.21, 0.01,  0.04,    0.06,  0.10]
para_s4['n'] =     [1.4,  1.6,   2.4,   3.4,    2.9,    3.4,  2.6,  1.4,   1.7,   2.1,   4.0,  3.9,   4.2,  1.8,   2.7,  2.4,   2.9,  2.6,     2.8,   2.8]

para_s4.to_csv("./para_s4.csv", index=False)
