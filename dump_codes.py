# dump
# # some test
# a, b, c, d = map(exprvar, 'abcd')
# a, b, c = map(bddvar, 'abc')
# f = a & b | a & c | b & c
# temp=f.to_dot()
#
# # Visualize the dot map
# s = Source(temp).view()
#
#
#
# def partition(input, nodes = 2):
#     # print("The first in is :",input[0])
#     # Nor input
#     out1 = Nor(input[0],input[1])
#     out2 = Nor(input[1],input[2])
#     out3 = Nor(input[0],input[2])
#
#
#     return
#
# def assembly(cell, QS =2):
#     return QS
#
# assembly(2,)
#
#
# expr("~a & b | ~c & d")
#
# # Give 3 input boolean Variable a,b,c
# a, b, c = map(exprvar, 'abc')
# # final = Nor(~b,Nor(a,b) ).to_dnf()
#
# # Cell with only 2 nodes ( nor , not )
# partition([a,b,c])
#
#
#
#
#
#
# e,f,g = map(bddvar, 'efg')
# Nor(a,b,c, simplify=True)
#
#
#
#
# # Indexing Variables
# x_0 = exprvar('x', 0)
# x_1 = exprvar('x', 1)
# x_0, x_1
# X = exprvars('x', 8)
# X
# # multi-dimensional
# X = exprvars('x', 4, 4)
# X
# a, b, c, d, e = map(exprvar, "abcde")
# f16 = Nor(a, b, c)
# f16.cardinality
# X = exp
# rvars('x', 4)
# f1 = truthtable(X, "0000011111010101")
# f1
#
# espresso_tts(f1)
# list(f1.satisfy_all())
# list(iter_points([a, b, c, d, e]))
# a, b = map(bddvar, 'ab')
# f = ~a & ~b | ~a & b | a & ~b | a & b ^a
# f
#
#
# f1 = truthtable(X, "0000011111010101")
# f2 = truthtable(X, "0000011111011101")
#
# a, b, c, d, e = map(exprvar, "abcde")
# Nor(a,b)
