import pandas as pd
import json
from pandas.io.json import json_normalize
import itertools
import numpy as np
from pyeda.inter import *
from graphviz import Source
import re
import pyeda.boolalg.expr as expr
import logic_fun as lf
import matplotlib.pyplot as plt
import time

# cd Circuit_bu # This line for ATOM-hydrogen path

# ---------------------------  Importing BESS 4 inputs NAND database -----------------------------
df4 = pd.read_csv("BESS_database/BESS_NAND.csv")
df4 = df4.replace(r'\r\n', ',', regex=True)

# --------- generating 4inputs database -------
Num_total_fun, Boolean_dict = lf.Boolean_espresso_gens(4)
convertion = dict([('Or', 'And'), ('And', 'Or')])



start = time.time()
"Waiting about 25 min"
Database = []
for func, expression in Boolean_dict.items():
    # converted expression
    func_converted = expr.expr(
        lf.multiple_replace(convertion, str(expression[0])))
    Num_in = len(func_converted.inputs)
    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gens(func_converted, 4, df4)
    data = np.append(np.append([func, expression], Opt_NOR_net), Num_in)
    Database.append(list(data))

end = time.time()
print("The Database generating time is {:.2f}s".format(end - start))


NOR4_database = np.array(Database)
NOR_df4 = pd.DataFrame({'Binary': NOR4_database[:, 0], 'Expresso_tts': np.array([i[0] for i in NOR4_database[:, 1]]), 'Num_gates': [
                       i[0] for i in NOR4_database[:, 2]], 'Optimal_net': [i[1] for i in NOR4_database[:, 2]], 'Permu_int': NOR4_database[:, 3], 'Num_in': NOR4_database[:, 4]})

# NOR_df4.to_csv("NOR_df4.csv", index=False)

set1 = NOR_df4[NOR_df4.Num_in == 1]
set2 = NOR_df4[NOR_df4.Num_in == 2]
set3 = NOR_df4[NOR_df4.Num_in == 3]
set4 = NOR_df4[NOR_df4.Num_in == 4]

set_list2 = set2.Expresso_tts.tolist()
set_list3 = set3.Expresso_tts.tolist()


expr_set2 = [expr.expr(i) for i in set_list2]
expr_set3 = [expr.expr(i) for i in set_list3]

expr_set23 = expr_set2+expr_set3
len(expr_set23)


df_orgin_combined = [df2, df3]
df_orgin_combined
t = 0
NOR_23_correct = []
for each in expr_set23:
    print('Number {} input expression is {}:'.format(t, each))
    var_list = list(each.inputs)
    vars_count = len(var_list)
    var_list_str = [str(i) for i in var_list]
    all_permu = list(itertools.permutations(var_list_str, vars_count))  # 💚
    origin = var_list_str
    pm_dict_list = [dict(zip(origin, i)) for i in all_permu]

    from math import factorial
    for per_num in range(factorial(vars_count)):  # 💚
        # set --------------------- each permutation ----------------
        each_permu = pm_dict_list[per_num]
        print('The permuations is:', each_permu)
        rs_ex_permu = expr.expr(lf.multiple_replace(each_permu, str(each)))
        rs_ex_permu_tt = list(rs_ex_permu.iter_relation())
        permu_br = ''.join([str(i[1]) for i in rs_ex_permu_tt])
        permu_intr = int(permu_br, 2)
        print('Int number:', permu_intr)

        df = df_orgin_combined[vars_count-2]  # 💚
        if not df[df['Integer'] == permu_intr].empty:

            opt_net = df[df['Integer'] == permu_intr]
            opt_struct = opt_net.values[0][3]
            opt_num_gates = opt_net.values[0][2]

            # generate reverse permuation mapping
            pm_reverse = dict((v, k) for k, v in pm_dict_list[per_num].items())
            print('reverse permutation is:', pm_reverse)
            # get reverse permuation in terms of abcd
            # convert_abcd = {'x[0]':'a','x[1]':'b','x[2]':'c','x[3]':'d'} # make the two variables to be a and b
            new_convert = dict(
                zip(each_permu.keys(), ['a', 'b', 'c'][0:vars_count]))  # 💚
            pm_reverse_new = dict((new_convert[key], new_convert[value]) for (
                key, value) in pm_reverse.items())
            print('reverse permuation final:', pm_reverse_new)

            NOR_net = lf.multiple_replace(
                pm_reverse_new, opt_struct).replace('NAND', 'NOR')
            print("original NAND design:", opt_struct)
            print("NOR design is:\n", NOR_net)

            NOR = [opt_num_gates, NOR_net]
            data = np.append([each, vars_count], NOR)
            NOR_23_correct.append(list(data))
            break
    t += 1
    # if t >1000: break
    print('\n')


# corrected set1 with 1 input
set1.Expresso_tts.str.contains('~')
set1['Num_gates'] = set1["Expresso_tts"].apply(lambda x: len(x)-4)

# for corrected df of 2,3 inputs
set23 = np.array(NOR_23_correct)
Binary = set2.append(set3).Binary
tts = [i[0] for i in set23]
Num_in = [i[1] for i in set23]
Num_gates = [i[2] for i in set23]
Optimal_net = [i[3] for i in set23]

# The correct df for input1, 2,3 and 4
dfc_1 = set1.drop(columns=['Permu_int'])
dfc_23 = pd.DataFrame({'Binary': Binary, 'Expresso_tts': tts,
                       'Num_gates': Num_gates, 'Optimal_net': Optimal_net, 'Num_in': Num_in})
dfc_4 = set4.drop(columns=['Permu_int'])


# The final corrected NOR database
dfc_final = dfc_1.append(dfc_23).append(dfc_4)

dfc_final.to_csv("NOR_database_Corrected_withoutID.csv", index=False)
