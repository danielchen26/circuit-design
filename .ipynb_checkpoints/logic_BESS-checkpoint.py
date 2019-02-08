import pandas as pd
import json
from pandas.io.json import json_normalize
import itertools
import numpy as np
from pyeda.inter import *
from graphviz import Source
import re, pyeda.boolalg.expr as expr
import logic_fun as lf


# ---------------------------  Importing BESS NAND database -----------------------------
# library of optimial gates
df = pd.read_csv("database/BESS_NAND.csv")
df = df.replace(r'\r\n', ',', regex = True)# modify df[['Optimal Network']] column


# ------------------- Generating a dictionary of minimal dnf design of each "Num_fun" of boolean functions
Num_fun, Boolean_dict = lf.Boolean_expresso_gen()
# convertion dict which switches 'Or' and 'And'
convertion = dict([('Or', 'And'),('And', 'Or')])

t = 1
NOR_database = []
for func, expression in Boolean_dict.items():
    print("The {}\n boolean function {}\n has minimal tts: {}".format(t, func, expression) )
    
    # converted expression
    func_converted = expr.expr(lf.multiple_replace(convertion, str(expression[0])))
    

    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gen(func_converted,df)

    
    data = np.append(func,Opt_NOR_net)
    NOR_database.append(data)
    t+=1
    if t >10:
        break




# # test example
# rs = str(Boolean_dict['1110101011111111'][0])
# Boolean_dict['1110101011111111'][0]
#
#
# # converted expression
# convertion = dict([('Or', 'And'),('And', 'Or')])
# rs_converted = expr.expr(lf.multiple_replace(convertion, rs))
#
# # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
# lf.permu_search_gen(rs_converted,df)
#
#
#
#
# # result truthtable (not used)
# rs_tt = expr2truthtable(rs_result, df)
#
# # truthtable list
# tt_list = list (rs_result.iter_relation())
# # truthtable binary Represenation
# bbr = ''.join([str(i[1]) for i in tt_list])
# # truthtable integer Represenation
# intr = int(bbr,2)
# intr
#
#
# # this example find a optimal network
# opt_net = df[df['Integer'] == intr]
# opt_net['Optimal Network']
#
#
# # -----------------------  Converting optimal Network expression in term of NOR gate
# # Corresponds to the specific permuation we convert the optimal Network back using that permutation
# opt_struct = opt_net.values[0][3]
# opt_struct.replace('NAND','NOR')
