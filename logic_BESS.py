import pandas as pd
import json
from pandas.io.json import json_normalize
import itertools
import numpy as np
from pyeda.inter import *
from graphviz import Source
import re, pyeda.boolalg.expr as expr
import logic_fun as lf


# ---------------------------  From BESS database -----------------------------
# library of optimial gates
df = pd.read_csv("~/Downloads/convertcsv.csv")
df = df.replace(r'\r\n', ',', regex = True)
df[['Optimal Network']]



# --------------------------   Creating all boolean functions --------------------------------
# -------------creating all possible binary output --------------
# ---output binary---
lst = list(itertools.product([0, 1], repeat=16))
n_lst = ''.join([str(i) for i in lst[0]])

test_tt = []
for each_tp in lst:
    new_tp = ''.join([str(i) for i in each_tp])
    test_tt.append(new_tp)


out_b = np.array(test_tt)
len(out_b)# total number possible function

# ---input binary ----
# in_b = ttvars('x', 4)
a,b,c,d= map(exprvar, "abcd")
in_b = farray([a,b,c,d])


# Logic minimization is known to be an NP-complete problem. It is equivalent to ﬁnding a minimal-cost set of subsets of a set 𝑆 that covers 𝑆. This is sometimes called the “paving problem”, because it is conceptually similar to ﬁnding the cheapest conﬁguration of tiles that cover a ﬂoor. Due to the complexity of this operation, PyEDA uses a C extension to the famous Berkeley Espresso library 1 .

# -- espresso_tts function to ﬁnd a low-cost, equivalent Boolean expression.
Boolean_dict = dict()
for each_outb in out_b:
    fi = truthtable(in_b, each_outb)
    fim = espresso_tts(fi)
    Boolean_dict[each_outb] = fim
Boolean_dict
# Data = {'Binary': Boolean_dict.keys(), 'Boolean_fun': Boolean_dict.values()}
# pd.DataFrame.from_dict(Data, orient='index')



ex1 = Boolean_dict['1110101011111111'][0].to_ast()
int('1110101011111111',2)





exchange = dict([('Or', 'And'),
                     ('And', 'Or')])

# test example
rs = str(Boolean_dict['1110101011111111'][0])
Boolean_dict['1110101011111111'][0]




# converted expression
rs_result = expr.expr(lf.multiple_replace(exchange, rs))



# Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
lf.permu_search_gen(rs_result,df)




# result truthtable (not used)
rs_tt = expr2truthtable(rs_result, df)

# truthtable list
tt_list = list (rs_result.iter_relation())
# truthtable binary Represenation
bbr = ''.join([str(i[1]) for i in tt_list])
# truthtable integer Represenation
intr = int(bbr,2)
intr


# this example find a optimal network
opt_net = df[df['Integer'] == intr]
opt_net['Optimal Network']


# -----------------------  Converting optimal Network expression in term of NOR gate
# Corresponds to the specific permuation we convert the optimal Network back using that permutation
opt_struct = opt_net.values[0][3]
opt_struct.replace('NAND','NOR')
