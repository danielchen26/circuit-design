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

# cd Circuit Scaling up(bottom-up) # This line for ATOM-hydrogen path

# ---------------------------  Importing BESS 4 inputs NAND database -----------------------------
df4 = pd.read_csv("BESS_database/BESS_NAND.csv")
df4 = df4.replace(r'\r\n', ',', regex=True)

# --------- generating 4inputs database -------
Num_total_fun, Boolean_dict = lf.Boolean_espresso_gens(4)
convertion = dict([('Or', 'And'), ('And', 'Or')])
Database = []


import time
start = time.time()
"the code you want to test stays here"
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
