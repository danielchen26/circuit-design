#Author: Tianchi Chen
# using DifferentialEquations, ParameterizedFunctions,Plots;gr()
using OrdinaryDiffEq
using ModelingToolkit
using Latexify,Random,Base
using CSV, DataFrames, Images
using Plots; pyplot()
##

# Define type
# mutable struct SType{T} <: DEDataVector{T}
#     x::Array{T,1}
#     p::T
# end

# mutable struct Time
#     p:: Float64
#     wait:: Float64
# end

mutable struct Param_E
    repressor
    RBS
    min:: Float64
    max:: Float64
    K:: Float64
    n:: Float64
end

# Import gate Parameters
para = CSV.read("DEmodels/param_db/para_s4.csv");
HlyIIR = Param_E(Matrix(para)[7,:]...)
LmrA = Param_E(Matrix(para)[10,:]...)
SrpR = Param_E(Matrix(para)[20,:]...)
BM3R1 = Param_E(Matrix(para)[6,:]...)
PhIF = Param_E(Matrix(para)[11,:]...)
PsrA = Param_E(Matrix(para)[14,:]...)
BetI = Param_E(Matrix(para)[3,:]...)


## -------  view response curve
# response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
# function response_crv(x, gate::Param_E)
#     min = gate.min; max = gate.max; K = gate.K; n = gate.n;
#     y = response.(min,max,K,n,x)
#     p = plot(x, y, xlims =(0,6),ylims = (0,1))
#     return p
# end
#
# response_crv(collect(1:0.1:6),HlyIIR)
# response_crv(collect(1:0.1:6),LmrA)
# response_crv(collect(1:0.1:6),SrpR)
# response_crv(collect(1:0.1:6),BM3R1)
# response_crv(collect(1:0.1:6),PhIF)
# response_crv(collect(1:0.1:6),PsrA)
# response_crv(collect(1:0.1:6),BetI)
# # -------  view response curve
γ = 0.025; ξ = 0.025
# response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
# degradation(x) = γ*x
# SR = @ode_def SR_latch begin
#     dm_HlyIIR = ξ*response(HlyIIR.min, HlyIIR.max, HlyIIR.K, HlyIIR.n, m_BetI + p) - degradation(m_HlyIIR)
#     dm_LmrA = ξ*response(LmrA.min, LmrA.max, LmrA.K, LmrA.n, m_HlyIIR + p) - degradation(m_LmrA)
#     dm_SrpR = ξ*response(SrpR.min, SrpR.max, SrpR.K, SrpR.n, m_HlyIIR + m_BetI) - degradation(m_SrpR)
#     dm_BM3R1 = ξ*response(BM3R1.min, BM3R1.max, BM3R1.K, BM3R1.n, m_LmrA + m_SrpR) - degradation(m_BM3R1)
#     dm_PhIF = ξ*response(PhIF.min, PhIF.max, PhIF.K, PhIF.n, m_BM3R1) - degradation(m_PhIF)
#     dm_PsrA = ξ*response(PsrA.min, PsrA.max, PsrA.K, PsrA.n, m_PhIF + m_BetI ) - degradation(m_PsrA)
#     dm_BetI = ξ*response(BetI.min, BetI.max, BetI.K, BetI.n, m_BM3R1 + m_PsrA) - degradation(m_BetI)
# end
#  ------- Define with MK
hill(dn, up, K, n, x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
### Define a differential equation system
@parameters t γ ξ p
@variables m_HlyIIR(t) m_LmrA(t) m_SrpR(t) m_BM3R1(t) m_PhIF(t) m_PsrA(t) m_BetI(t)
@derivatives D'~t
eqs1 = [
    # Bit 1 =================
    D(m_HlyIIR) ~ ξ * hill(HlyIIR.min, HlyIIR.max, HlyIIR.K, HlyIIR.n, m_BetI + p)        - deg(m_HlyIIR),
    D(m_LmrA ) ~ ξ * hill(LmrA.min, LmrA.max, LmrA.K, LmrA.n, m_HlyIIR + p)       - deg(m_LmrA),
    D(m_SrpR  ) ~ ξ * hill(SrpR.min, SrpR.max, SrpR.K, SrpR.n, m_HlyIIR + m_BetI) - deg(m_SrpR),
    D(m_BM3R1 ) ~ ξ * hill(BM3R1.min, BM3R1.max, BM3R1.K, BM3R1.n, m_LmrA + m_SrpR)   - deg(m_BM3R1),
    D(m_PhIF) ~ ξ * hill(PhIF.min, PhIF.max, PhIF.K, PhIF.n, m_BM3R1)            - deg(m_PhIF),
    D(m_PsrA ) ~ ξ * hill(PsrA.min, PsrA.max, PsrA.K, PsrA.n, m_PhIF + m_BetI) - deg(m_PsrA),
    D(m_BetI ) ~ ξ * hill(BetI.min, BetI.max, BetI.K, BetI.n, m_BM3R1 + m_PsrA)  - deg(m_BetI)]
de1 = ODESystem(eqs1, t, [m_HlyIIR, m_LmrA, m_SrpR, m_BM3R1, m_PhIF, m_PsrA,m_BetI], [γ,ξ,p])
counter_ecoli = ODEFunction(de1)

## Figure 3: Eq dynamics
u0 = rand(1:22., length(counter_ecoli.syms))
p = 0.0
tspan = (0.0, 1000.0)
param = [0.025,0.025,p]
prob0 = ODEProblem(counter_ecoli, u0, tspan, param)
sol0 = solve(prob0, Tsit5())
plot(sol0, lw = 2, ylims = (0, 22),
     xlabel = " Time steps", ylabel = "Concentration",
     color = [:blue :red],
     vars=[:m_PsrA,:m_BetI],
     label =["PsrA" "BetI"],
     title = "Randomized initial conditions yield stable states")

plt = plot()
for i = 1:50
    u0 = rand(1:22., length(counter_ecoli.syms))
    prob0 = ODEProblem(counter_ecoli, u0, tspan, param)
    sol0 = solve(prob0, Tsit5())
    plot!(plt, sol0, lw = 1, ylims = (0, 22),
         xlabel = " Time steps", ylabel = "Concentration",
         vars=[:m_PsrA,:m_BetI],
         color = [:blue :orange],
         # label =["PsrA" "BetI"],
         label = "",
         legend = false,
         title = "Randomized initial conditions yield stable states",
         dpi = 300)
end
plt
sol0 = solve(prob0, Tsit5())
plot!(sol0,vars=[:m_PsrA,:m_BetI],
      color = [:blue :orange],
      label =["Gate 6: PsrA" "Gate 7: BetI"],
      xlabel = " Time steps", ylabel = "Concentration",
      legend = true,la=0.6)
savefig(plt,"two steady states.png")
## --`




function run_prob_1bit(;init_relax, duration,relax,signal,K,n,up)
    u0 =  rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = 20)
    param = [up,0.002,K,n,0.025,0.025,p]
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










































##
# Random.seed!(134)
# u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
# p = 0.0
# prob0 = ODEProblem(SR,u0,(0.0,1000.0),p)
# sol0 = solve(prob0,Tsit5())
# plot(sol0,vars=[:m_HlyIIR,:m_BetI,:m_PsrA])


# u0_1 = SType(sol0[end].x, 20.0)
# p=20.0
# prob1 = ODEProblem(SR,u0_1,(0.0,2200.0),p)
# sol1 = solve(prob1,Tsit5())
# plot(sol1,vars=[(0,6),(0,7)], lw =2,xlabel = "time", ylabel = "concentration")



# # using cbs -------------- testing ---------
#
# SR = @ode_def_bare SR_latch begin
#     dm_HlyIIR = ξ*response(HlyIIR.min, HlyIIR.max, HlyIIR.K, HlyIIR.n, m_BetI + p) - degradation(m_HlyIIR)
#     dm_LmrA = ξ*response(LmrA.min, LmrA.max, LmrA.K, LmrA.n, m_HlyIIR + p) - degradation(m_LmrA)
#     dm_SrpR = ξ*response(SrpR.min, SrpR.max, SrpR.K, SrpR.n, m_HlyIIR + m_BetI) - degradation(m_SrpR)
#     dm_BM3R1 = ξ*response(BM3R1.min, BM3R1.max, BM3R1.K, BM3R1.n, m_LmrA + m_SrpR) - degradation(m_BM3R1)
#     dm_PhIF = ξ*response(PhIF.min, PhIF.max, PhIF.K, PhIF.n, m_BM3R1) - degradation(m_PhIF)
#     dm_PsrA = ξ*response(PsrA.min, PsrA.max, PsrA.K, PsrA.n, m_PhIF + m_BetI ) - degradation(m_PsrA)
#     dm_BetI = ξ*response(BetI.min, BetI.max, BetI.K, BetI.n, m_BM3R1 + m_PsrA) - degradation(m_BetI)
# end p
#
# u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
# p = [0.0]
# tspan = (0.0,5000.0)
# ts, cb = make_cb2([1000,3200],1,20.)
# prob = ODEProblem(SR,u0,tspan,p)
# sol = solve(prob,Tsit5(),callback=cb, tstops=ts)
#
#
# function make_cb2(ts_in, var1, var1_val)
#     ts = ts_in
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#       if integrator.t == ts[1]
#           integrator.p[var1] = var1_val
#       elseif integrator.t == ts[2]
#           integrator.p[var1] = 0.0
#       end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     return ts, cb
# end








# Find the lb and ub of the P inpluse and the min time to wait for the next signal

# Find the local minimum index
locs2 =  findlocalminima(sol1[6,:])
ind_min1 = [i[1] for i in locs2][2]

#  Test the time starting from ti =1 to ti = first local min
t_wait = []
ul_range =[]
t_range = []
anim = @animate for ti = 1: ind_min1
    u0_t11 = sol1[:,ti]
    println("The duration of P is:", sol1.t[ti])
    p=0
    prob11 = ODEProblem(SR,u0_t11,(0.0,3000.0),p)
    sol11 = solve(prob11,Tsit5())
    display(plot(sol11,vars=[(0,6),(0,7)]))
#      How much time to wait until the system become statble again
    locs =findlocalmaxima(sol11[6,:])[1]
    display(plot!(sol11[6,:]))
    stable_t_ind = locs[1]
    if sol11[7,end] < sol11[6,end]
        push!(ul_range, ti)
        push!(t_range, sol1.t[ti])
        push!(t_wait, sol11.t[stable_t_ind])
    end
end
gif(anim,fps = 2)

#  give the min time to wait(lower bound) for the next signal within the switching range
maximum(t_wait)
# switching index range
ul_range
# switching lower bound time
t_lb = t_range[1]
# switching upper bound time
t_ub = t_range[end]



# For loop simulating more switching cycles

# set P duration= 600, and t_wait = 1500

# Initialize vector of the final plot
P_set = [];sol_set =[];t_set = []

#  1. make system goes to steady state
rng = MersenneTwister(124)
p=0; u0= Float64[i for i in rand(rng,1:22,7)]
prob_steady = ODEProblem(SR,u0,(0.0,1500.0),p); sol_steady = solve(prob_steady,SSRootfind())
plot(sol_steady,vars=[(0,6),(0,7)])

p1 = zeros(size(sol_steady.t))
push!(t_set, sol_steady.t)

#  Set the time parameters
time = Time(400, 570)

for cycle = 1:1
#   1. Add inpulse p=20 for 600
    p = 20;
    if cycle ==1
        u0_1= SType(sol_steady[end], 20.0)
    else
        u0_1= SType(sol_relax[end].x, 20.0)
    end
    prob_impulse = ODEProblem(SR,u0_1,(0.0,time.p),p) ; sol_impulse = solve(prob_impulse,Tsit5())
    display(plot(sol_impulse,vars=[(0,6),(0,7)]))

    p_on = zeros(size(sol_impulse.t)) .+ 20
    if cycle ==1
        t1_n = sol_impulse.t .+ sol_steady.t[end]
        push!(t_set, t1_n)
    else
        t1_n = sol_impulse.t .+ t2_n[end]
        push!(t_set, t1_n)
    end

#   2. Relax system for 1500
    p = 0; u0_2= SType(sol_impulse[end].x, 20.0)
    prob_relax = ODEProblem(SR,u0_2,(0.0,time.wait),p) ; sol_relax = solve(prob_relax,Tsit5())
    display(plot(sol_relax,vars=[(0,6),(0,7)]))

    p_off = zeros(size(sol_relax.t))
    t2_n = sol_relax.t .+ t1_n[end]
    push!(t_set, t2_n);

    push!(P_set, p_on); push!(P_set, p_off);
    push!(sol_set, [i.x for i in sol_impulse.u]);   push!(sol_set, [i.x for i in sol_relax.u]);
end

t_set_f = vcat(t_set...)
p_set_f = vcat([p1; P_set]...)
sol_set_f = [sol_steady.u;vcat(sol_set...)]

PsrA = [i[6] for i in sol_set_f]
BetI = [i[7] for i in sol_set_f]

#  Final plot
plot(t_set_f, PsrA,lw =2)
plot!(t_set_f, BetI, marker = (:star,2),lw = 2)
plot!(t_set_f, p_set_f, lw = 2, line = (:dot, :arrow, 1.9, 3, :green), xlabel = "time", ylabel = "concentration")
title!("Simulation Result")
