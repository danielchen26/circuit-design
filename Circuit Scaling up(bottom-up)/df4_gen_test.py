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

# ---------------------------  Importing BESS 4 inputs NAND database -----------------------------
df4 = pd.read_csv("BESS_database/BESS_NAND.csv")
df4 = df4.replace(r'\r\n', ',', regex=True)


# ------------------- Generating a dictionary of minimal dnf design of each "Num_fun" of boolean functions
Num_fun, Boolean_dict = lf.Boolean_espresso_gen()

# convertion dict which switches 'Or' and 'And'
convertion = dict([('Or', 'And'), ('And', 'Or')])

# Creating NOR optimal net database(takes quite a bit of time!!!)
database = []
t = 1
for func, expression in Boolean_dict.items():
    print("The {}\n boolean function {}\n has minimal tts: {}".format(
        t, func, expression))

    # converted expression
    func_converted = expr.expr(
        lf.multiple_replace(convertion, str(expression[0])))

    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gen(func_converted, df4)

    data = np.append([func, expression], Opt_NOR_net)
    database.append(list(data))
    t += 1
    if t > 10:
        break


NOR_database = np.array(database)

NOR_df = pd.DataFrame({'Binary': NOR_database[:, 0], 'Expresso_tts': np.array(
    [i[0] for i in NOR_database[:, 1]]), 'Num_gates': NOR_database[:, 2], 'Optimal_net': NOR_database[:, 3]})
