# randomly assign simulated gates for 1-bit counter.
#  Modified from silico_sim.jl 1bit section

#Author: Tianchi Chen
## ------ Import package and functions
using ModelingToolkit, OrdinaryDiffEq #DifferentialEquations
using Plots; gr(fontfamily = "Souce Code Pro for Powerline"); #pyplot()#
using Latexify, Random, Base
using CSV, DataFrames, ProgressMeter
include("functions.jl")# using BlackBoxOptim, LinearAlgebra
##

## ==== Build multiple counter connectors : 1Bit counter case 📗 =========
# ==== Define ODEProblem =======
# γ = 0.025
# ξ = 0.025
# up = 1.5;
# dn = 0.002;
# K = 0.081;
# n = 2.81;
mutable struct Hill{T}
	dn::T
    up::T
    K::T
    n::T
end

# @unpack dn, up, K, n = Hill(0.002, 1.5, 0.081, 2.81)
mono_param = Hill(0.002, 1.5, 0.081, 2.81)


# hill(param::Hill, x) = param.dn + (param.up - param.dn) * param.K^(param.n) / (param.K^(param.n) + abs(x)^(param.n))
hill(dn, up, K, n, x) = dn + (up - dn) * K^(n) / (K^(n) + abs(x)^(n))
deg(x) = γ * x
#-------------------------- Define a differential equation system
# @parameters t up dn K n γ ξ p
@parameters t γ ξ p
@parameters LexA1_dn LexA1_up LexA1_K LexA1_n 
@parameters IcaR_dn IcaR_up IcaR_K IcaR_n
@parameters CI1_dn CI1_up CI1_K CI1_n
@parameters PsrA_dn PsrA_up PsrA_K PsrA_n
@parameters BM3RI_dn BM3RI_up BM3RI_K BM3RI_n
@parameters HKCI_dn HKCI_up HKCI_K HKCI_n
@parameters PhlF_dn PhlF_up PhlF_K PhlF_n
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
@derivatives D'~t
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(LexA1_dn, LexA1_up, LexA1_K, LexA1_n, m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(IcaR_dn, IcaR_up, IcaR_K, IcaR_n, m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(CI1_dn, CI1_up, CI1_K, CI1_n, m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(PsrA_dn, PsrA_up, PsrA_K, PsrA_n, m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(BM3RI_dn, BM3RI_up, BM3RI_K, BM3RI_n, m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(HKCI_dn, HKCI_up, HKCI_K, HKCI_n, m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(PhlF_dn, PhlF_up, PhlF_K, PhlF_n, m1_PsrA + m1_HKCI)  - deg(m1_PhlF)]
de1 = ODESystem(eqs1, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF], 
[ LexA1_dn, LexA1_up, LexA1_K, LexA1_n, 
  IcaR_dn, IcaR_up, IcaR_K, IcaR_n,
  CI1_dn, CI1_up, CI1_K, CI1_n, 
  PsrA_dn, PsrA_up, PsrA_K, PsrA_n, 
  BM3RI_dn, BM3RI_up, BM3RI_K, BM3RI_n, 
  HKCI_dn, HKCI_up, HKCI_K, HKCI_n, 
  PhlF_dn, PhlF_up, PhlF_K, PhlF_n,
  γ, ξ, p])
ode_f1 = ODEFunction(de1)

## Run the 1bit counter problem
# # ==== Randomized Initials equlibrations =====
u0 = rand(1:22., length(ode_f1.syms))
p = 0.0
tspan = (0.0, 3000.0)
mono_set = [0.002, 1.5, 0.081, 2.81]
param = [ mono_set...,
          0.002, 1.5, 0.081, 2.81,
          0.002, 1.5, 0.081, 2.81,
          0.002, 1.5, 0.081, 2.81,
          0.002, 1.5, 0.081, 2.81,
          0.002, 1.5, 0.081, 2.81,
          0.002, 1.5, 0.081, 2.81, 
          0.025,0.025,p]
prob0 = ODEProblem(ode_f1, u0, tspan, param)
sol0 = solve(prob0, Tsit5())
plot(sol0, lw = 2, ylims = (0, 22))











mutable struct gate_param_assign{T}
	g1::T
    g2::T
    g3::T
    g4::T
    g5::T
    g6::T
    g7::T
end

## Run one example
function run_prob_1bit(;init_relax, duration,relax,signal, gate_p_set::gate_param_assign)
    u0 =  rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = 20)
    param = [ gate_p_set.g1...,
              gate_p_set.g2...,
              gate_p_set.g3...,
              gate_p_set.g4...,
              gate_p_set.g5...,
              gate_p_set.g6...,
              gate_p_set.g7...,
              0.025,0.025,p]
    prob0 = ODEProblem(ode_f1, u0, tspan, param)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
    return sol, ts
end
# #  ==== Induced by signals with t_on time ======
# # when constant signal is applied, we expected to see oscillation
# p = 20.0;
# t_on = 1000.0;
# prob1 = ODEProblem(silico, sol0[end], (0.0, t_on), p)
# sol1 = solve(prob1, SSRootfind())
# plot(sol1, vars = [:m_HKCI, :m_PhlF], lw = 2)

sol, ts = run_prob_1bit(;init_relax = 2500., duration=270.,relax=2500.,
                        signal=20.,K=0.081,n=2.81, up= 1.5);
up= 1.5
py = plot(sol, vars = [:m1_HKCI,:m1_PhlF],
          lw = 1.5, xlabel = "time", ylabel = "concentration",
          title = "Signal Duration: 270", ylims = (0.,3*up))
##
plt = plot( sol,vars = [:m1_HKCI,:m1_PhlF],lw = 1.5,
            xlabel = "time", ylabel = "concentration",
            title = "Signal Duration: 270",
            ylims = (0.,3*up))
scatter!(plt,ts,ones(length(ts)))

# plot(sol,vars = [:m1_HKCI,:m1_PhlF], tspan=(40000.0,60000.0)) # plot the dynamics of at time [40000.0,60000.0]
## ======== Sampling Parameters for 1 bit counter ===================
using ProgressMeter
# # Varying K, n, δ
# df = DataFrame(K = Float64[], n = Float64[], δ =  Float64[])
# tt = []
# @time @showprogress for K = 0.001: 0.01:0.1, n = 1.5:0.1:3., δ = 250:5:350
#   K = K; n = n;
#   sol = run_prob(duration = δ)
#   # display(plot(sol, vars = [:m_HKCI, :m_PhlF]))
#   # Check switching =============
#   g6i_av, g6ii_av, g7i_av, g7ii_av = switching(sol,4)
#
#   if g6i_av > g7i_av
#       dn < g6ii_av < (up+dn)/2  &&  (up+dn)/2 < g7ii_av < up ? push!(df, [K, n, δ]) : nothing
#   elseif g6i_av < g7i_av
#       dn < g7ii_av < (up+dn)/2  &&   (up+dn)/2 < g6ii_av < up ? push!(df, [K, n, δ]) : nothing
#   end
#
#   push!(tt,1.)
#   @show sum(tt)
#   # if sum(tt) >= 200.
#   #     break
#   # end
# end
# df

# function Carrying_bit(sol, ts)
#     ts_ID = []
#     C1 = []; C2 = [];
#     for ti = 1:2:length(ts)
#         push!(ts_ID, ti, ti+1)
#         # connector 1
#         opt = L_max(sol,10, ts[ti],ts[ti + 1])
#         lmax1 = -opt.minimum;
#
#         # connector 2
#         opt2 = L_max(sol,20, ts[ti],ts[ti + 1])
#         lmax2 = -opt2.minimum;
#
#         push!(C1,lmax1)
#         push!(C2,lmax2)
#     end
#     return ts_ID, C1, C2
# end
# ts_ID, C1, C2  = Carrying_bit(sol, ts)

# function MB_T0(ts_ID, C1, C2; thred =up/2)
#     # Identify the T0
#     # C1_B = C1 .> thred; C2_B = C2 .> thred
#     T0 = []
#     ts_cycle = collect(1:2:length(ts_ID))
#     for i = 1:length(C1)
#         if C1[i] > thred && C2[i] > thred
#             # @show i , C1[i], C2[i], thred
#             push!(T0,  ts_cycle[i], ts_cycle[i] + 1)
#             break
#         else
#             continue
#         end
#     end
#
#     return T0
# end
# # Identify the T0 for multibit counter
# T0 = MB_T0(ts_ID, C1, C2; thred =up/2) # if something goes wrong, this may give me 0

function Switch(sol, idx::Array , Δ0, ts, t_id::Array )
    # Get the ID for gates from input array
    ID1 = idx[1]; ID2 = idx[2];

    # Show eq sol for two output
    G6 = sol(Δ0)[ID1];  G7 = sol(Δ0)[ID2]
    # @show G6, G7
    g6i_tot = g7i_tot = 0

    T_afs = (ts[t_id[2]] + ts[t_id[2] + 1])/2
    G6_afs = sol(T_afs)[ID1]; G7_afs = sol(T_afs)[ID2]

    return G6_afs, G7_afs
end

function Switch_cost(sol, idx::Array, Δ0, ts, T0)
    # G6 is HKCI, G7 is PhlF
    G6 =[];G7 =[];
    for sw in 0:2:14
        B1_1, B1_2 = Switch(sol,idx, Δ0, ts, T0 .+ sw)
        # println("Bit: ", B1_1, B1_2)
        push!(G6, B1_1);     push!(G7, B1_2)
    end
    return G6, G7
end

function osci(G6,up)
    thred = up/2
    G6b = G6 .> thred
    G6b1 = [G6b[i] for i = 1:2:length(G6b)]
    G6b2 = [G6b[i] for i = 2:2:length(G6b)]
    @show G6b1, G6b2, G6b1 .+ G6b2
    return sum(G6b1 .+ G6b2), sum(G6b1 .* G6b2) # comp should be 4, and sum should be 0
end


function cost_bit1(sol, ts, up)
    G6, G7 = Switch_cost(sol, [6,7],1000.,ts, [1,2])
    comp6, dot6 = osci(G6,up); comp7, dot7 = osci(G7,up);
    if comp6 ==4 && comp7 == 4 && dot6 == 0 && dot7 == 0
        println("good")
        return costtot = 0
    else
        println("bad")
        return costtot = 1
    end
end

cost_bit1(sol, ts, up)


# run the searching
df_1bit = DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[], up = Float64[])
tt = []
@time @showprogress for K = 0.01: 0.05:1., n = 1.:0.2:10., δ = 10:5:500, A = 20., up = 1:1.:10 #@showprogress
    # for K = 0.011: 0.05:1., n = 1.5:0.5:10., δ = 250:10:350, A = 20., up = 1:0.5:10 # try broader range
    # K = 0.011: 0.02:0.1, n = 1.5:0.1:3., δ = 250:5:350, A = 20., up = 1:0.1:3, # original
    # solve DiffEq given parameters
    sol, ts = run_prob_1bit(;init_relax = 2500., duration=δ, relax=2000., signal=A, K=K, n=n, up = up)
    # Get cost for this parameters set
    costtot = cost_bit1(sol, ts, up)
    println("K:$K n:$n, δ:$δ, A:$A, up:$up\n")
    @show costtot
    costtot == 0 ? push!(df_1bit, [K, n, δ, A, up]) : nothing

    # count example
    push!(tt,1.)
    @show sum(tt)
    # if sum(tt) >= 800.
    #     break
    # end
end

CSV.write("1Bit_DB_(K = 0.01: 0.05:1., n = 1.:0.2:10., δ = 10:5:500, A = 20., up = 1:1.:10).csv", df_1bit)


# import 1bit database
using CSV, DataFrames
df_1b = CSV.read("1Bit_DB.csv")
df_1b[10000,:].K, df_1b[10000,:].n, df_1b[10000,:].up


Sample(4, 1:150)
gate_p_set = gate_param_assign([0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n],
                               [0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n],
                               [0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n],
                               [0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n],
                               [0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n],
                               [0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n],
                               [0.002, df_1b[10000,:].up, df_1b[10000,:].K, df_1b[10000,:].n])

sol, ts = run_prob_1bit(;init_relax = 2500., duration=320.,relax=2500., signal=20., gate_p_set);
plot(sol)

