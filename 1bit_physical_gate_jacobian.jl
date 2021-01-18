##
using ModelingToolkit,OrdinaryDiffEq
using CSV, DataFrames
mutable struct gate_param_assign{T}
	g1::T
    g2::T
    g3::T
    g4::T
    g5::T
    g6::T
    g7::T
end

db = CSV.read("1Bit_DB_(K = 0.01: 0.05:1., n = 1.:0.2:10., δ = 10:5:500, A = 20., up = 1:1.:10).csv",DataFrame)
# load para_s4 table
df = CSV.read("./DEmodels/param_db/para_s4.csv",DataFrame)
rename!(df, [:repressor, :RBS, :dn, :up, :K, :n], makeunique=true)

##
f(dn, up, K, n, x) = dn + (up - dn) * K^(n) / (K^(n) + abs(x)^(n))
g(x) = γ * x

@parameters t γ ξ p
@parameters LexA1[1:4] IcaR[1:4] CI1[1:4] PsrA[1:4] BM3RI[1:4] HKCI[1:4] PhlF[1:4]
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
@derivatives D'~t
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * f(LexA1..., m1_PhlF + p)            - g(m1_LexA1),
    D(m1_IcaR ) ~ ξ * f(IcaR...,  m1_LexA1 + p)           - g(m1_IcaR),
    D(m1_CI1  ) ~ ξ * f(CI1...,   m1_LexA1 + m1_PhlF)     - g(m1_CI1),
    D(m1_PsrA ) ~ ξ * f(PsrA...,  m1_IcaR + m1_CI1)       - g(m1_PsrA),
    D(m1_BM3RI) ~ ξ * f(BM3RI..., m1_PsrA)                - g(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * f(HKCI...,  m1_BM3RI + m1_PhlF)     - g(m1_HKCI),
    D(m1_PhlF ) ~ ξ * f(PhlF...,  m1_PsrA + m1_HKCI)      - g(m1_PhlF)]
de1 = ODESystem(eqs1, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF],
[ LexA1..., IcaR..., CI1..., PsrA..., BM3RI..., HKCI..., PhlF..., γ, ξ, p])
ode_f1 = ODEFunction(de1)

## solve equations
p_unique = Any[[0.07, 4.3, 0.05, 1.7], [0.2, 5.9, 0.19, 1.8], [0.03, 2.8, 0.21, 2.4], [0.2, 2.2, 0.18, 2.1], [0.07, 2.5, 0.19, 2.6], [0.07, 3.8, 0.41, 2.4], [0.06, 3.8, 0.07, 1.6]]
# the cello_set
p_unique = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
gate_p_set = gate_param_assign(p_unique...)
# test equilibrium for the example
u0 =  rand(1:22., length(ode_f1.syms))
p = 0;tspan = (0.0, 3000.0)
param = [ gate_p_set.g1...,
       gate_p_set.g2...,
       gate_p_set.g3...,
       gate_p_set.g4...,
       gate_p_set.g5...,
       gate_p_set.g6...,
       gate_p_set.g7...,
       0.025,0.025,p]
prob0 = ODEProblem(ode_f1, u0, tspan, param)
sol0 = solve(prob0, Tsit5() )
plot(sol0)

## calculate jacobian in two ways?🔴🔴
using ForwardDiff
jm = ForwardDiff.jacobian(u -> ode_f1(u,param,tspan[end]),sol0[end])
using LinearAlgebra
eigen(jm)

# Question 2:🔺
# how to use generate_jacobian and build_function to build jacobian function that can be evaluated with any given u and p?
# build_function(eval(generate_jacobian(de1)[2]),[m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF],
# [ LexA1..., IcaR..., CI1..., PsrA..., BM3RI..., HKCI..., PhlF..., γ, ξ, p])
jac = eval(generate_jacobian(de1)[2],[m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF],[LexA1..., IcaR..., CI1..., PsrA..., BM3RI..., HKCI..., PhlF..., γ, ξ, p])
jac