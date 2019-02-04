import pandas as pd
import re
import pyeda.boolalg.expr as expr
import itertools
import numpy as np
from pyeda.inter import *




# Switching 'Or' and 'And'
def multiple_replace(adict, text):
  import pyeda.boolalg.expr as expr
  # Create a regular expression from all of the dictionary keys
  regex = re.compile("|".join(map(re.escape, adict.keys(  ))))
  return regex.sub(lambda match: adict[match.group(0)], text)




# # this step we do a permutation and recorded in a dictionary
# Here is to write a function for doing all permutation and convert all permuation string back to pyeda expression
def permu_search_gen(bf,df):
    '''
    bf: input of the PyEDA binary function (simplified by espresso_tts);
    df: the external library loaded for the origin NAND gate 
    '''
    import pyeda.boolalg.expr as expr
    all_permu = list(itertools.permutations('abcd',4))
    origin = ('a','b','c','d')
    pm_dict_list = [dict(zip(origin, each)) for each in all_permu]

    # assinging input boolean function
    rs_result = bf

    # The following block test each permutation case in the permutation dictionary list,
    # totally 24 possible case, and choose one the match integer that corresponds to the data library
    for per_num in range(24):
        # test --------------------- each permutation ----------------
        #
        rs_exchange = multiple_replace(pm_dict_list[per_num], str(rs_result))

        # 3 possible correction terms, only one is needed to be corrected.
        correction = {'Anc': 'And', 'Anb': 'And', 'Ana': 'And'}
        for key in correction.keys():
            if key in rs_exchange:
                print(key)
                correction_sub = {key: correction[key]}
                rs_ex_permu = expr.expr(multiple_replace(correction_sub,rs_exchange))
            else:
                rs_ex_permu = expr.expr(multiple_replace(correction,rs_exchange))


        # rs_ex_permu =  expr.expr(multiple_replace(correction,rs_exchange))
        rs_ex_permu_tt = list (rs_ex_permu.iter_relation())
        permu_br = ''.join([str(i[1]) for i in rs_ex_permu_tt])
        permu_intr = int(permu_br,2)
        # permu_intr
        if not df[df['Integer'] == permu_intr].empty:
            opt_net = df[df['Integer'] == permu_intr]
            opt_struct = opt_net.values[0][3]
            print('The original network structure is :\n', opt_struct)
            print('\nThe lib origin term has been found:\n')
            print('we have use the permuation dictionary of ', pm_dict_list[per_num])

            # generate reverse permuation mapping
            pm_reverse = dict((v,k) for k, v in pm_dict_list[per_num].items())
            print('The reversed permutation mapping is:\n', pm_reverse)

            # construct the new NOR gate library
            NOR_net = multiple_replace(pm_reverse, opt_struct).replace('NAND','NOR')
            print("The optimal NOR network structure is:\n", NOR_net)

            break# break the searching permuation loops of 24 possible cases
        else:
            print('This specific permu is not in the lib')
        #
        # test --------------------- each permutation ----------------
