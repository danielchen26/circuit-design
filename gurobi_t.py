# gurobi
from gurobipy import *
import collections
import numpy as np
# tested with Python 2.7.6 & Gurobi 6.5
maxi = range(1,7+1)
maxj = range(1,2+1)
maxk = range(1,4+1)
logcr = range(4,7+1)
alpha = collections.defaultdict(int, {(2,2): 1, (3,1): 1, (4,1): 1, (4,2): 1}) # will return 0 if key is not in dict
innor_indexes = {1: (2,3), 2: (4,5), 3: (6,7)}
innor_limit = 2
logb_indexes = [(1,2), (1,3), (2,4), (2,5), (3,6), (3,7)]
output = collections.defaultdict(int, {2:1, 3:1})
output_gate = 1


model = Model('Logic Design')
model.Params.UpdateMode = 1
in1 = {}
out = {}
nor = {}
for i in maxi:
    for j in maxj:
        in1[i,j] = model.addVar(vtype=GRB.BINARY, name="in1_i{}_j{}".format(i,j))
    for k in maxk:
        out[i,k] = model.addVar(vtype=GRB.BINARY, name="out_i{}_k{}".format(i,k))
    nor[i] = model.addVar(vtype=GRB.BINARY, name="nor_i{}".format(i))


# Constraints
# enor
for i in maxi:
    for j in maxj:
        model.addConstr(in1[i,j] <= nor[i], name="enor_i{}_j{}".format(i,j))


# innor
for i, (g1, g2) in innor_indexes.iteritems():
    model.addConstr(nor[g1] + nor[g2] + quicksum(in1[i,j] for j in maxj) <= innor_limit, name="innor_i{}".format(i))







# loga
for k in maxk:
    for (g1, g2) in logb_indexes:
        model.addConstr(out[g1,k] + out[g2,k] <= 1, name="loga_k{}_g{}_g{}".format(k,g1,g2))
# lc
for k in maxk:
    for i in maxi:
        for j in maxj:
            model.addConstr(alpha[k,j]*in1[i,j]+out[i,k] <= 1, name="lc_i{}_j{}_k{}".format(i,j,k))
# logc
for k in maxk:
    for i, (g1, g2) in innor_indexes.iteritems():
        model.addConstr(out[g1,k] + out[g2,k] + quicksum(alpha[k,j]*in1[i,j] for j in maxj) + out[i,k] - nor[i] >= 0,
                        name="logc_k{}_i{}".format(k,i))
    for i in logcr:
        model.addConstr(quicksum(alpha[k,j]*in1[i,j] for j in maxj) + out[i,k] - nor[i] >= 0,
                        name="logc_k{}_i{}".format(k,i))
