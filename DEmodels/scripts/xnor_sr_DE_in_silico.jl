#Author: Tianchi Chen
using DifferentialEquations, ParameterizedFunctions,Plots;gr()
using Latexify,Random,Base
using CSV, DataFrames, Images
using Parameters
include("functions.jl")
# Define type
mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::T
end

mutable struct Time
    p:: Float64
    wait:: Float64
end

@with_kw mutable struct Hill
    up:: Float64
    dn:: Float64
    K:: Float64
    n:: Float64
end


γ = 0.025
ξ = 0.025
hill(P::Hill,x) = (P.dn + (P.up - P.dn)*P.K^P.n/(P.K^P.n + x^P.n))
degradation(x) = γ*x

P = Hill(up = 1.5, dn =0.002, K= 0.11, n =2.81)

silico = @ode_def_bare counter begin
    dm_LexA1 = ξ*hill(P, m_PhlF + p)        - degradation(m_LexA1)
    dm_IcaR  = ξ*hill(P, m_LexA1 + p)       - degradation(m_IcaR)
    dm_CI1   = ξ*hill(P, m_LexA1 + m_PhlF)  - degradation(m_CI1)
    dm_PsrA  = ξ*hill(P, m_IcaR + m_CI1)    - degradation(m_PsrA)
    dm_BM3RI = ξ*hill(P, m_PsrA)            - degradation(m_BM3RI)
    dm_HKCI  = ξ*hill(P, m_BM3RI + m_PhlF ) - degradation(m_HKCI)
    dm_PhlF  = ξ*hill(P, m_PsrA + m_HKCI)   - degradation(m_PhlF)
end p


u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
p = 0.0
prob0 = ODEProblem(silico,u0,(0.0,1000.0),p)
sol0 = solve(prob0,SSRootfind())
plot(sol0, lw = 2, ylims = (0,22))
# plot!(fig, sol0, vars =[:m_HKCI,:m_PhlF], label = ["HKCI" "PhlF"], lw = 2, ylims = (0,22), linecolor = [:orange :green], legend = true)




u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
p = 0.0
prob0 = ODEProblem(silico,u0,(0.0,10000.0),p)
sol0 = solve(prob0,SSRootfind())
py_ss = plot(sol0,vars=[:m_HKCI,:m_PhlF],linecolor = [:orange :green])#

u0_1 = SType(sol0[end].x, 20.0)
p=20.0
prob1 = ODEProblem(silico,u0_1,(0.0,5200.0),p)
sol1 = solve(prob1,SSRootfind())
py = plot(sol1,vars=[:m_HKCI,:m_PhlF], lw =2,xlabel = "time", ylabel = "concentration")
# savefig(py,"~/Desktop/silico_oscilation.png")




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
    prob11 = ODEProblem(silico,u0_t11,(0.0,3000.0),p)
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









# =========== For loop simulating more switching cycles =================
# set P duration= 600, and t_wait = 1500

# Initialize vector of the final plot
P_set = []; sol_set =[]; t_set = []

#  1. make system goes to steady state
rng = MersenneTwister(124)
p=[0]; u0= Float64[i for i in rand(rng,1:22,7)]
prob_steady = ODEProblem(silico,u0,(0.0,1500.0),p); sol_steady = solve(prob_steady,SSRootfind())
plot(sol_steady,vars=[:m_HKCI,:m_PhlF])

p1 = zeros(size(sol_steady.t))
push!(t_set, sol_steady.t)

#  Set the time parameters
time = Time(458.7, 1500)


for cycle = 1:2
#   1. Add inpulse p=20 for 600
    p = [20.];
    if cycle ==1
        u0_1= SType(sol_steady[end], p[1])
    else
        u0_1= SType(sol_relax[end].x, p[1])
    end
    prob_impulse = ODEProblem(silico,u0_1,(0.0,time.p),p);
    sol_impulse = solve(prob_impulse,SSRootfind())
    # display(plot(sol_impulse,vars=[(0,6),(0,7)]))

    p_on = zeros(size(sol_impulse.t)) .+ p[1]
    if cycle ==1
        t1_n = sol_impulse.t .+ sol_steady.t[end]
        push!(t_set, t1_n)
    else
        t1_n = sol_impulse.t .+ t2_n[end]
        push!(t_set, t1_n)
    end

#   2. Relax system for 1500
    p = [0]; u0_2= SType(sol_impulse[end].x, 0.0)
    prob_relax = ODEProblem(silico,u0_2,(0.0,time.wait),p);
    global sol_relax = solve(prob_relax,SSRootfind())
    # display(plot(sol_relax,vars=[(0,6),(0,7)]))

    p_off = zeros(size(sol_relax.t))
    global t2_n = sol_relax.t .+ t1_n[end]
    push!(t_set, t2_n)

    push!(P_set, p_on);     push!(P_set, p_off);
    push!(sol_set, [i.x for i in sol_impulse.u]);      push!(sol_set, [i.x for i in sol_relax.u]);
end

t_set_f = vcat(t_set...)
p_set_f = vcat([p1; P_set]...)
sol_set_f = [sol_steady.u;vcat(sol_set...)]

PsrA = [i[6] for i in sol_set_f]
BetI = [i[7] for i in sol_set_f]

#  Final plot
plot(t_set_f, PsrA,lw =2, label ="HKCI")
plot!(t_set_f, BetI, marker = (:star,2),lw = 2, label ="PhlF")
plot!(t_set_f, p_set_f, lw = 2, line = (:dot, :arrow, 1.9, 3, :green), xlabel = "time", ylabel = "concentration", label = "Signal")
title!("Simulation Result")












# =============== Wraping above functions =============

function action_bound()
    u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
    p = 0.0
    prob0 = ODEProblem(silico,u0,(0.0,10000.0),p)
    sol0 = solve(prob0,SSRootfind())
    # py_ss = plot(sol0,vars=[:m_HKCI,:m_PhlF],linecolor = [:orange :green])#

    u0_1 = SType(sol0[end].x, 20.0)
    p=20.0
    prob1 = ODEProblem(silico,u0_1,(0.0,5200.0),p)
    sol1 = solve(prob1,SSRootfind())
    # py = plot(sol1,vars=[:m_HKCI,:m_PhlF], lw =2,xlabel = "time", ylabel = "concentration")
    # savefig(py,"~/Desktop/silico_oscilation.png")


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
        # println("The duration of P is:", sol1.t[ti])
        p=0
        prob11 = ODEProblem(silico,u0_t11,(0.0,3000.0),p)
        sol11 = solve(prob11,Tsit5())
        # display(plot(sol11,vars=[(0,6),(0,7)]))
    #      How much time to wait until the system become statble again
        locs =findlocalmaxima(sol11[6,:])[1]
        # display(plot(sol11[6,:]))
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
    return t_lb,t_ub
end




t_lb,t_ub


action_bound()







using ProgressMeter
Σ = 0
@time @showprogress for up = 0.5:0.5:1.5, dn = 0.:0.0005:0.002, K = 0.01:0.04:0.11, n = 2.:0.458:2.81
    @show up, dn, K, n
    # redefine Models
    P = Hill(up = up, dn =dn, K= K, n =n)
    silico = @ode_def_bare counter begin
        dm_LexA1 = ξ*hill(P, m_PhlF + p)        - degradation(m_LexA1)
        dm_IcaR  = ξ*hill(P, m_LexA1 + p)       - degradation(m_IcaR)
        dm_CI1   = ξ*hill(P, m_LexA1 + m_PhlF)  - degradation(m_CI1)
        dm_PsrA  = ξ*hill(P, m_IcaR + m_CI1)    - degradation(m_PsrA)
        dm_BM3RI = ξ*hill(P, m_PsrA)            - degradation(m_BM3RI)
        dm_HKCI  = ξ*hill(P, m_BM3RI + m_PhlF ) - degradation(m_HKCI)
        dm_PhlF  = ξ*hill(P, m_PsrA + m_HKCI)   - degradation(m_PhlF)
    end p

    # Initialize vector of the final plot
    P_set = []; sol_set =[]; t_set = []

    #  1. make system goes to steady state
    rng = MersenneTwister(124)
    p=[0]; u0= Float64[i for i in rand(rng,1:5,7)]
    prob_steady = ODEProblem(silico,u0,(0.0,1500.0),p); sol_steady = solve(prob_steady,SSRootfind())
    plot(sol_steady,vars=[:m_HKCI,:m_PhlF])

    p1 = zeros(size(sol_steady.t))
    push!(t_set, sol_steady.t)

    #  Set the time parameters

    time = Time(458.7, 1500)
    for cycle = 1:4
    #   1. Add inpulse p=20 for 600
        p = [5.];
        if cycle ==1
            u0_1= SType(sol_steady[end], p[1])
        else
            u0_1= SType(sol_relax[end].x, p[1])
        end
        prob_impulse = ODEProblem(silico,u0_1,(0.0,time.p),p);
        sol_impulse = solve(prob_impulse,SSRootfind())
        # display(plot(sol_impulse,vars=[(0,6),(0,7)]))

        p_on = zeros(size(sol_impulse.t)) .+ p[1]
        if cycle ==1
            t1_n = sol_impulse.t .+ sol_steady.t[end]
            push!(t_set, t1_n)
        else
            t1_n = sol_impulse.t .+ t2_n[end]
            push!(t_set, t1_n)
        end

    #   2. Relax system for 1500
        p = [0]; u0_2= SType(sol_impulse[end].x, 0.0)
        prob_relax = ODEProblem(silico,u0_2,(0.0,time.wait),p);
        global sol_relax = solve(prob_relax,SSRootfind())
        # display(plot(sol_relax,vars=[(0,6),(0,7)]))

        p_off = zeros(size(sol_relax.t))
        global t2_n = sol_relax.t .+ t1_n[end]
        push!(t_set, t2_n)

        push!(P_set, p_on);     push!(P_set, p_off);
        push!(sol_set, [i.x for i in sol_impulse.u]);      push!(sol_set, [i.x for i in sol_relax.u]);
    end

    t_set_f = vcat(t_set...)
    p_set_f = vcat([p1; P_set]...)
    sol_set_f = [sol_steady.u;vcat(sol_set...)]

    PsrA = [i[6] for i in sol_set_f]
    BetI = [i[7] for i in sol_set_f]

    #  Final plot
    pp = plot()
    plot!(pp,t_set_f, PsrA,lw =2, label ="HKCI")
    plot!(pp,t_set_f, BetI, lw = 2, label ="PhlF") # marker = (:star,2),
    plot!(pp,t_set_f, p_set_f, lw = 2,  xlabel = "time", ylabel = "concentration", label = "Signal") # line = (:dot, :arrow, 1.9, 3, :green),
    title!("Simulation Result \n Max: $up  Min: $dn  K: $K   n:$n")
    display(pp)
    global Σ = Σ + 1
    println("Cycle:", Σ)
end
