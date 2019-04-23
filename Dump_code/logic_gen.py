import numpy as np
import pandas as pd
import scipy as sp
from pyeda.inter import *
from graphviz import Source
import itertools

# --------------------------   Main --------------------------------
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
possible = 0
for each_outb in out_b:
    fi = truthtable(in_b, each_outb)
    fim = espresso_tts(fi)

    if fim[0].to_ast()[0] == 'or' and 'lit' not in [i[0] for i in fim[0].to_ast()]:
        ori = fim[0].to_ast()[1::]
        ary  =  list(map(list, ori))

        mat = []
        for i in range(len(ary)):
            mat.append([j[1] for j in ary[i][1::]])

        count = 0
        tt = np.array(mat)
        length = [len(i) for i in tt]
        tt_array = [np.array(i) for i in mat]
        p_len = [sum(i>0) for i in tt_array]

        cases = tuple(zip(length,p_len))
        if (4,4) in cases or (4,3) in cases:
            count += 1
        # print('next function\n')
        if count !=0:
            possible+=1
possible




possible2 = 0
for each_outb in out_b:
    fi = truthtable(in_b, each_outb)
    fim = espresso_tts(fi)

    if fim[0].to_ast()[0] != 'or' and len(fim[0].to_ast()[1::]) != 1:

        ori = fim[0].to_ast()[1::]
        mat2 = [y for x,y in ori]
        # mat2
        # ary
        # ary  =  list(map(list, ori))

        # mat2 = [each[1] for each in ary]
        # mat2
        if (len(mat2), sum(np.array(mat2)>0)) in [(4,4),(4,3)]:
            # print(mat2)
            possible2+=1
possible2


for each_outb in out_b[1:1000]:
    fi = truthtable(in_b, each_outb)
    fim = espresso_tts(fi)

    if 'lit' in [i[0] for i in fim[0].to_ast()] and len(fim[0].to_ast()[1::]) != 1:
        print("ok", fim)
        print([i[0] for i in fim[0].to_ast()].count('lit'))
    # else:
    #     print('not',fim)
