
#Author: Tianchi Chen
## ------ Import package and functions
using DifferentialEquations, ModelingToolkit
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
using Latexify, Random, Base
using CSV, DataFrames, ProgressMeter
include("functions.jl")# using BlackBoxOptim, LinearAlgebra
## ==== Build multiple counter connectors : 1Bit counter case 📗 =========
# ==== Define ODEProblem =======
# γ = 0.025
# ξ = 0.025
# up = 1.5;
# dn = 0.002;
# K = 0.081;
# n = 2.81;
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
 # Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
@derivatives D'~t
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF)]
de1 = ODESystem(eqs1, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f1 = ODEFunction(de1)

## Run the 1bit counter problem
# # ==== Randomized Initials equlibrations =====
u0 = rand(1:22., length(ode_f1.syms))
p = 0.0
tspan = (0.0, 3000.0)
param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
prob0 = ODEProblem(ode_f1, u0, tspan, param)
sol0 = solve(prob0, Tsit5())
plot(sol0, lw = 2, ylims = (0, 22))
# #
## Run one example
function run_prob_1bit(;init_relax, duration,relax,signal,K,n,up,cycle)
    u0 =  rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = cycle)
    param = [up,0.002,K,n,0.025,0.025,p]
    # prob0 = ODEProblem(ode_f1, u0, (0,init_relax), param)
    # sol0 = solve(prob0, Tsit5())
    prob = ODEProblem(ode_f1, u0, tspan, param)
    sol = solve(prob, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
end
# #  ==== Induced by signals with t_on time ======
# # when constant signal is applied, we expected to see oscillation
# p = 20.0;
# t_on = 1000.0;
# prob1 = ODEProblem(silico, sol0[end], (0.0, t_on), p)
# sol1 = solve(prob1, SSRootfind())
# plot(sol1, vars = [:m_HKCI, :m_PhlF], lw = 2)

sol, ts = run_prob_1bit(;init_relax = 5000., duration=270.,relax=5000.,signal=20.,K=0.081,n=2.81, up= 1.5,cycle=20);
up= 1.5
py = plot(sol, vars = [:m1_HKCI,:m1_PhlF],
          lw = 1.5, xlabel = "time", ylabel = "concentration",
          title = "Signal Duration: 270", ylims = (0.,3*up))
##
plt = plot(sol,vars = [:m1_HKCI,:m1_PhlF],lw = 1.5, xlabel = "time", ylabel = "concentration",title = "Signal Duration: 270",ylims = (0.,3*up))
scatter!(plt,ts,ones(length(ts)))

## ======== Sampling Parameters for 1 bit counter ===================
# using ProgressMeter

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
## check 1bit database
db1 = CSV.read("1Bit_DB.csv")
db1
for i =1#:size(db1)[1]
    # i = rand(1:size(db12)[1])
    δδ = db1[i,:].δ; AA =db1[i,:].A;  KK = db1[i,:].K; nn = db1[i,:].n; upp = db1[i,:].up
    sol, ts = run_prob_1bit(;init_relax = 5000., duration=δδ,relax=5000.,signal=AA, K=KK, n=nn, up= upp, cycle=20)
    param = [KK, nn, δδ, AA, upp]
    plt = plot(sol,vars = [:m1_HKCI,:m1_PhlF],lw = 1.5, xlabel = "time", ylabel = "concentration",title = "$param",ylims = (0.,3*up))
    scatter!(plt,ts,ones(length(ts)))
    cost_bit1(sol, ts, upp)
end

i =1

sol, ts = run_prob_1bit(;init_relax = 5000., duration=δδ,relax=5000.,signal=AA, K=KK, n=nn, up= upp, cycle=20)
param = [KK, nn, δδ, AA, upp]
plt = plot(sol,vars = [:m1_HKCI,:m1_PhlF],lw = 1.5, xlabel = "time", ylabel = "concentration",title = "$param",ylims = (0.,3*up))
scatter!(plt,ts,ones(length(ts)))
