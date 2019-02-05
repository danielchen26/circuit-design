# using DifferentialEquations, Plots
# using Latexify
using JuMP, Gurobi, Clp


# example

m = Model(solver = GurobiSolver())
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )

@objective(m, Max, 5x + 3*y )
@constraint(m, 1x + 5y <= 3.0 )

print(m)

status = solve(m)

println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))


#            Model Formulation
# -----------------------------------------

# I be a set of gates,
# J={A,B} be a set of two Inputs
# K={1,…,4} be the set of Input/Output-Indexes.


# We have given the limit of NOR inputs l

# -------Defining parameters
# For each gate i ∈ I¹, we have given the two input gates:
#   gᵢ¹  : first input gate
#   gᵢ²  : second input gate


# For each input j ∈ J and each Input/Output-Index k ∈ K, we have the parameter:
#   αₖⱼ  : value of the jth input in the kth row of the truth table


# For each Input/Output-Index k ∈ K, we have the parameter:
#   Oₖ  : output value of kth row of the truth table.


# ------ Defining Variables
# 1. each gate i has binary variable :
#        norᵢ ∈ {1, 0}   →    1 if ith gate exists.

# 2. each gate i ∈ I and each input j ∈ I , define the binary variable:
#        inᵢⱼ  : 1 if the jth input is an input to the ith gate.

# 3. For each gate i ∈ I and each Input/Output-Index k ∈ K, define the binary variable:
#        Oᵢₖ   : will take the value of the kth entry in the truth table.


# ------ Defining constraints
# 1. a NOR gate can only have external input if it exists:
#       inᵢⱼ ≤ norᵢ   ∀i ∈ I, j ∈ J

# 2. if a NOR gate has 1 or 2 external inputs leading into it, only 1 or no NOR gates can feed into it.
#       nor_gᵢ¹ + nor_gᵢ² + Σⱼ αₖⱼ * inᵢⱼ ≤ l  ∀i ∈ I¹

# 3. This constraints make sure, that the output signal from NOR gate i must be the correct logical function (NOR) of the input signals into gate i if it exists. The α describes for the jth row and for the kth external input the value of the external input.
#       O_gᵢ¹k + O_gᵢ²k ≤ 1     ∀ i ∈ I¹, k ∈ K
#       αₖⱼ * inᵢⱼ + Oᵢₖ ≤ 1      ∀ i ∈ I, k ∈ K
#
