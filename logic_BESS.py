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


# Creating NOR optimal net database(takes quite a bit of time!!!)

database = []
t = 1
for func, expression in Boolean_dict.items():
    print("The {}\n boolean function {}\n has minimal tts: {}".format(t, func, expression) )

    # converted expression
    func_converted = expr.expr(lf.multiple_replace(convertion, str(expression[0])))

    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gen(func_converted,df)

    data = np.append([func,expression],Opt_NOR_net)
    database.append(list(data))
    t+=1
    if t>10:break


NOR_database = np.array(database)


NOR_df = pd.DataFrame({'Binary':NOR_database[:,0],'Expresso_tts':np.array([i[0] for i in NOR_database[:,1]]),'Num_gates':NOR_database[:,2],'Optimal_net':NOR_database[:,3]})
NOR_df.to_csv("NOR_df.csv", index=False)




# ---------------------------------------------------------------------

# Test specific boolean function
st = '1111111100000000'
st = '0000001101010111'
express = Boolean_dict[st]# ~d
func_converted = expr.expr(lf.multiple_replace(convertion, str(express[0])))

import pyeda.boolalg.expr as expr
all_permu = list(itertools.permutations('abcd',4))
origin = ('a','b','c','d')
pm_dict_list = [dict(zip(origin, each)) for each in all_permu]

# assinging input boolean function
rs_result = func_converted

# The following block test each permutation case in the permutation dictionary list,
# totally 24 possible case, and choose one the match integer that corresponds to the data library

NOR = []
for per_num in range(24):
    #
    rs_exchange = lf.multiple_replace(pm_dict_list[per_num], str(rs_result))
    rs_exchange

    # 3 possible correction terms, only one is needed to be corrected.
    correction = {'Anc': 'And', 'Anb': 'And', 'Ana': 'And'}
    for key in correction.keys():
        if key in rs_exchange:
            print(key)
            correction_sub = {key: correction[key]}
            rs_ex_permu = expr.expr(lf.multiple_replace(correction_sub,rs_exchange))
        else:
            rs_ex_permu = expr.expr(lf.multiple_replace(correction,rs_exchange))


    # rs_ex_permu =  expr.expr(multiple_replace(correction,rs_exchange))
    rs_ex_permu_tt = list (rs_ex_permu.iter_relation())
    permu_br = ''.join([str(i[1]) for i in rs_ex_permu_tt])
    permu_intr = int(permu_br,2)

    if not df[df['Integer'] == permu_intr].empty:
            opt_net = df[df['Integer'] == permu_intr]
            opt_struct = opt_net.values[0][3]
            opt_num_gates = opt_net.values[0][2]
            print('The original network structure is :\n', opt_struct)

            print('we have use the permuation dictionary of ', pm_dict_list[per_num])

            # generate reverse permuation mapping
            pm_reverse = dict((v,k) for k, v in pm_dict_list[per_num].items())
            print('The reversed permutation mapping is:\n', pm_reverse)

            # construct the new NOR gate library
            NOR_net = lf.multiple_replace(pm_reverse, opt_struct).replace('NAND','NOR')
            print("The optimal NOR network structure is:\n", NOR_net)

            NOR = [opt_num_gates, NOR_net]
            break

    print(df[df['Integer'] == permu_intr])
# ---------------------------------------------------------------------



# ------- Some visualizations ------
import matplotlib.pyplot as plt
import altair as alt

NOR_reload =  pd.read_csv("NOR_df2.csv")
type(NOR_reload.Binary[3])
NOR_reload['Binary'] = list(Boolean_dict.keys())# somehow need to force the string type
columns = NOR_reload.columns

NOR_stat = NOR_reload.groupby('Num_gates').count().reset_index().replace('NAN',np.nan).dropna()
NOR_stat['Num_gates'] = NOR_stat['Num_gates'].astype(int)
NOR_stat.sort_values(by=['Num_gates'])#.apply(lambda x: x.cumsum())
# NOR_reload.to_csv("NOR_df2.csv", index =False)
NOR_stat
NOR_reload[NOR_reload['Binary'] =='1111111100000000']
