#Author: Tianchi Chen
using DifferentialEquations, ParameterizedFunctions,Plots;gr()
using Latexify,Random,Base
using CSV, DataFrames, Images
# Define type
mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::T
end

mutable struct Time
    p:: Float64
    wait:: Float64
end

mutable struct Param_Y
    max:: Float64
    min:: Float64
    K:: Float64
    n:: Float64
end

include("functions.jl")

# Import gate Parameters
para = CSV.read("./DEmodels/param_db/para_y.csv");# load yeast parameters

LexA1 = Param_Y(Matrix(para)[8,2:end]...)
IcaR = Param_Y(Matrix(para)[7,2:end]...)
CI1 = Param_Y(Matrix(para)[2,2:end]...)
PsrA = Param_Y(Matrix(para)[11,2:end]...)
BM3RI = Param_Y(Matrix(para)[1,2:end]...)
HKCI = Param_Y(Matrix(para)[6,2:end]...)
PhlF = Param_Y(Matrix(para)[10,2:end]...)

# # -------  view response curve
# response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
# function response_crv(x, gate::Param_Y)
#     min = gate.min; max = gate.max; K = gate.K; n = gate.n;
#     y = response.(min,max,K,n,x)
#     p = plot(x, y)
#     return p
# end
#
# st = 0; intv = 0.0001; fn = 3
# response_crv(collect(st:intv:fn),LexA1)
# response_crv(collect(st:intv:fn),IcaR)
# response_crv(collect(st:intv:fn),CI1)
# response_crv(collect(st:intv:fn),PsrA)
# response_crv(collect(st:intv:fn),BM3RI)
# response_crv(collect(st:intv:fn),HKCI)
# response_crv(collect(st:intv:fn),PhlF)
# # -------  view response curve

γ = 0.025
ξ = 0.025
response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
degradation(x) = γ*x



yeast = @ode_def_bare counter begin
    dm_LexA1 = ξ*response(LexA1.min, LexA1.max, LexA1.K, LexA1.n, m_PhlF + p) - degradation(m_LexA1)
    dm_IcaR = ξ*response(IcaR.min, IcaR.max, IcaR.K, IcaR.n, m_LexA1 + p) - degradation(m_IcaR)
    dm_CI1 = ξ*response(CI1.min, CI1.max, CI1.K, CI1.n, m_LexA1 + m_PhlF) - degradation(m_CI1)
    dm_PsrA = ξ*response(PsrA.min, PsrA.max, PsrA.K, PsrA.n, m_IcaR + m_CI1) - degradation(m_PsrA)
    dm_BM3RI = ξ*response(BM3RI.min, BM3RI.max, BM3RI.K, BM3RI.n, m_PsrA) - degradation(m_BM3RI)
    dm_HKCI = ξ*response(HKCI.min, HKCI.max, HKCI.K, HKCI.n, m_BM3RI + m_PhlF ) - degradation(m_HKCI)
    dm_PhlF = ξ*response(PhlF.min, PhlF.max, PhlF.K, PhlF.n, m_PsrA + m_HKCI) - degradation(m_PhlF)
end p

p = 0.0
# Random.seed!(134)

# ------ randomize Initials
fig = plot()
anim = @animate for rd = 1:8
    u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
    p = 0.0
    prob0 = ODEProblem(yeast,u0,(0.0,1000.0),p)
    sol0 = solve(prob0,SSRootfind())
    plot!(fig, sol0, vars =[:m_HKCI,:m_PhlF], lw = 2, ylims = (0,22), linecolor = [:orange :green],label =["HKCI" "PhlF"] )
end
gif(anim, "/tmp/tmp.gif", fps = 10)


u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
p = 0.0
prob0 = ODEProblem(yeast,u0,(0.0,10000.0),p)
sol0 = solve(prob0,Tsit5())
py_ss = plot(sol0,vars=[:m_HKCI,:m_PhlF],linecolor = [:orange :green])#

u0_1 = SType(sol0[end].x, 20.0)
p=20.0
prob1 = ODEProblem(yeast,u0_1,(0.0,2200.0),p)
sol1 = solve(prob1,Tsit5())
py = plot(sol1,vars=[:m_HKCI,:m_PhlF], lw =2,xlabel = "time", ylabel = "concentration")
# savefig(py,"~/Desktop/yeast_oscilation.png")







# # test with cbs  ---------------------------------
yeast_cb = @ode_def_bare counter begin
    dm_LexA1 = ξ*response(LexA1.min, LexA1.max, LexA1.K, LexA1.n, m_PhlF + p) - degradation(m_LexA1)
    dm_IcaR = ξ*response(IcaR.min, IcaR.max, IcaR.K, IcaR.n, m_LexA1 + p) - degradation(m_IcaR)
    dm_CI1 = ξ*response(CI1.min, CI1.max, CI1.K, CI1.n, m_LexA1 + m_PhlF) - degradation(m_CI1)
    dm_PsrA = ξ*response(PsrA.min, PsrA.max, PsrA.K, PsrA.n, m_IcaR + m_CI1) - degradation(m_PsrA)
    dm_BM3RI = ξ*response(BM3RI.min, BM3RI.max, BM3RI.K, BM3RI.n, m_PsrA) - degradation(m_BM3RI)
    dm_HKCI = ξ*response(HKCI.min, HKCI.max, HKCI.K, HKCI.n, m_BM3RI + m_PhlF ) - degradation(m_HKCI)
    dm_PhlF = ξ*response(PhlF.min, PhlF.max, PhlF.K, PhlF.n, m_PsrA + m_HKCI) - degradation(m_PhlF)
end p


u0 = SType(Float64[i for i in rand(1:22.,7)], 0.0)
p = [0.0]
tspan = (0.0,5000.0)
# ts, cb = make_cb([1000,1600],1,20.)
# ts, cb =make_cb2([1000,1600,3000,3600],  20., 0., 20., 0.)
ts, cb = mk_cb_mul([1000,1600,3000,3600], 20., 0., 20., 0.)
prob = ODEProblem(yeast_cb,u0,tspan,p)
sol = solve(prob,SSRootfind(),callback=cb, tstops=ts, reltol = 1e-15,abstol = 1e-19)
plot(sol,vars=[:m_HKCI,:m_PhlF])


using Interact
u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
p = [0.0]
@manipulate for Δ = 1.:1000., δ = 1.:50.
    ts, cb = mk_cb_mul([1000, 1000 + Δ, 3000 , 3000 + Δ], δ, 0., δ, 0.)
    prob = ODEProblem(yeast_cb,u0,tspan,p)
    sol = solve(prob,SSRootfind(),callback=cb, tstops=ts, reltol = 1e-20,abstol = 1e-30)
    # plot(sol)
    plot(sol,vars=[:m_HKCI,:m_PhlF])
end


# Find the lb and ub of the P inpluse and the min time to wait for the next signal
# Find the local minimum index
locs2 =  findlocalminima(sol1[6,:])
ind_min1 = [i[1] for i in locs2][2]

#  Test the time starting from ti =1 to ti = first local min
t_wait = []
ul_range =[]
t_range = []
for ti = 1: ind_min1
    u0_t11 = sol1[:,ti]
    println("The duration of P is:", sol1.t[ti])
    p=0
    prob11 = ODEProblem(yeast,u0_t11,(0.0,3000.0),p)
    sol11 = solve(prob11,Tsit5())
    display(plot(sol11,vars=[(0,6),(0,7)]))
#      How much time to wait until the system become statble again
    locs =findlocalmaxima(sol11[6,:])[1]
    display(plot(sol11[6,:]))
    stable_t_ind = locs[1]
    if sol11[7,end] < sol11[6,end]
        push!(ul_range, ti)
        push!(t_range, sol1.t[ti])
        push!(t_wait, sol11.t[stable_t_ind])
    end
end


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
P_set = []; sol_set =[]; t_set = [];

#  1. make system goes to steady state
rng = MersenneTwister(124)
p=[0]; u0= Float64[i for i in rand(rng,1:22,7)]
prob_steady = ODEProblem(yeast,u0,(0.0,1500.0),p); sol_steady = solve(prob_steady,SSRootfind())
plot(sol_steady,vars=[(0,6),(0,7)])

p1 = zeros(size(sol_steady.t))
push!(t_set, sol_steady.t)

#  Set the time parameters
time = Time(350, 1500)


for cycle = 1:4
#   1. Add inpulse p=20 for 600
    p = [20];
    if cycle ==1
        u0_1= SType(sol_steady[end], 20.0)
    else
        u0_1= SType(sol_relax[end].x, 20.0)
    end
    prob_impulse = ODEProblem(yeast,u0_1,(0.0,time.p),p) ; sol_impulse = solve(prob_impulse,Tsit5())
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
    p = [0]; u0_2= SType(sol_impulse[end].x, 0.0)
    prob_relax = ODEProblem(yeast,u0_2,(0.0,time.wait),p);
    global sol_relax = solve(prob_relax,Tsit5())
    display(plot(sol_relax,vars=[(0,6),(0,7)]))

    p_off = zeros(size(sol_relax.t))
    global t2_n = sol_relax.t .+ t1_n[end]
    push!(t_set, t2_n)

    push!(P_set, p_on);     push!(P_set, p_off);
    push!(sol_set, [i.x for i in sol_impulse.u]);      push!(sol_set, [i.x for i in sol_relax.u]);
end

t_set_f = vcat(t_set...)
p_set_f = vcat([p1; P_set]...)
sol_set_f = [sol_steady.u;vcat(sol_set...)]

HKCI_sol = [i[6] for i in sol_set_f];   PhlF_sol = [i[7] for i in sol_set_f]

#  Final plot
plot(t_set_f, HKCI_sol,lw =2)
plot!(t_set_f, PhlF_sol, marker = (:star,2),lw = 2)
plot!(t_set_f, p_set_f, lw = 2, line = (:dot, :arrow, 1.9, 3, :green), xlabel = "time", ylabel = "concentration")
title!("Simulation Result")
















# # =================== Cost function ===================
# Θ(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
