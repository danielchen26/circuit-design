
#Author: Tianchi Chen
## ------ Import package and functions
using ModelingToolkit, OrdinaryDiffEq #DifferentialEquations
using Plots; pyplot()#gr(fontfamily = "Souce Code Pro for Powerline");
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
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
### Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
@derivatives D'~t
##
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF)]
##
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


##


# # ===== Animation of control problems =======
# sig_rg = 250:10:330
# anim = Anim_gen_controlproblem(cycle, Δ0, Δ, δ, A, sig_rg)
# gif(anim, "/tmp/anim.gif", fps = 5)


##
# Find the lb and ub of the P inpluse and the min time to wait for the next signal
# Find the local minimum index
#  ========== Mark the local minimum ==========

# ====== Local minimum and maximum with t_on ============
# using Images
# t_on = 1000.0
# plt_min, mark6_min, mark7_min = localminmax_sig(sol, t_on, "lmin")
# plt_min
# mark7_min
# plt_max, mark6_max, mark7_max = localminmax_sig(sol, t_on, "lmax")
# plt_max







##  Test the time starting from ti =1 to ti = first local min
# # t_wait = []
# ul_range = []
# t_range = []
# for ti = 1:mark7_min[2]
#     u0_t11 = sol[:, ti]
#     println("Signal Duration : ", sol.t[ti])
#     p = 0
#     param = [1.5,0.002,0.81,2.81,0.025,0.025,p]
#     prob11 = ODEProblem(ode_f1, u0_t11, (0.0, 3000.0), param)
#     sol11 = solve(prob11, Tsit5())
#     display(plot(sol11, vars = [:m1_HKCI, :m1_PhlF]))
#
#     # How much time to wait until the system become statble again
#     locs = findlocalmaxima(sol11[7, :])#[1]
#     marklist = [i[1] for i in locs]
#     marker_spec = (:circle, 1, 0.6, :purple, stroke(2, 0.5, :orange, :dot))
#     display(scatter!(sol11[marklist], vars = [:m1_PhlF], xlims = (1:t_on), marker = marker_spec))
#     # display(plot(sol11[6,:]))
#     # stable_t_ind = locs[1]
#     if sol11[7, end] > sol11[6, end]
#         push!(ul_range, ti)
#         push!(t_range, sol.t[ti])
#         # push!(t_wait, sol11.t[stable_t_ind])
#     end
# end
#
# @show t_range
# @show ul_range
#
#
# ti = 45
# u0_t11 = sol1[:, ti]
# println("Signal Duration : ", sol1.t[ti])
# p = 0
# prob11 = ODEProblem(silico, u0_t11, (0.0, 3000.0), p)
# sol11 = solve(prob11, Tsit5())
# display(plot(sol11, vars = [(0, 6), (0, 7)]))
#
# locs = findlocalmaxima(sol11[7, :])#[1]
# marklist = [i[1] for i in locs]
# marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
# display(scatter!(sol11[marklist], vars = [:m_PhlF], xlims = (1:t_on), marker = marker_spec))
#
#
#
#
#
# #  give the min time to wait(lower bound) for the next signal within the switching range
# maximum(t_wait)
# # switching index range
# ul_range
# # switching lower bound time
# t_lb = t_range[1]
# # switching upper bound time
# t_ub = t_range[end]
## need to rewrite 🔴


# ------ bboptimize way of searching Parameters
# function cost(pa)
#     # try to optimize δ: signal duration, K: dissociation rate, n: hill para
#     up = pa[1]
#     K = pa[2]
#     n = pa[3]
#     @show Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control()# δ = param
#     prob0 = ODEProblem(silico, u0, tspan, p)
#     sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
#     norm(sol[7, end] - up)
# end
# bound = [(1.5, 5.5), (1.5, 5.5), (1.5, 5.5)]
# result = bboptimize(cost; SearchRange = bound, NumDimensions = 3)


# ========  To run a problem 1 time, given input of signal duration



























## ======= Build multiple counter connectors : 2 Bits counter case 📗 =========
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
### Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t) g1(t) g2(t) g3(t) m2_LexA1(t) m2_IcaR(t) m2_CI1(t) m2_PsrA(t) m2_BM3RI(t) m2_HKCI(t) m2_PhlF(t)
@derivatives D'~t
##
eqs2 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF),
    # Connector 1 =============
    D(g1)      ~ ξ * hill(p)                   - deg(g1),
    D(g2)      ~ ξ * hill(m1_HKCI)             - deg(g2),
    D(g3)      ~ ξ * hill(g1 + g2)             - deg(g3), # g3 sserves as the input for the 2nd bit
    # Bit 2 =============
    D(m2_LexA1) ~ ξ * hill(m2_PhlF + g3)       - deg(m2_LexA1),
    D(m2_IcaR ) ~ ξ * hill(m2_LexA1 + g3)      - deg(m2_IcaR),
    D(m2_CI1  ) ~ ξ * hill(m2_LexA1 + m2_PhlF) - deg(m2_CI1),
    D(m2_PsrA ) ~ ξ * hill(m2_IcaR + m2_CI1)   - deg(m2_PsrA),
    D(m2_BM3RI) ~ ξ * hill(m2_PsrA)            - deg(m2_BM3RI),
    D(m2_HKCI ) ~ ξ * hill(m2_BM3RI + m2_PhlF) - deg(m2_HKCI),
    D(m2_PhlF ) ~ ξ * hill(m2_PsrA + m2_HKCI)  - deg(m2_PhlF)]
##
de2 = ODESystem(eqs2, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF, g1, g2, g3, m2_LexA1, m2_IcaR, m2_CI1, m2_PsrA,m2_BM3RI, m2_HKCI, m2_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f2 = ODEFunction(de2)
# ##
# u0 = Float64[i for i in rand(1:22, 17)]
# Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 270, cycle = 10)
# prob0 = ODEProblem(silico2, u0, tspan, p)
# sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)

## Run the 2bit counter problem
# # ==== Randomized Initials equlibrations =====
u0 = rand(1:22., length(ode_f2.syms))
p = 0.0
tspan = (0.0, 3000.0)
param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
prob0 = ODEProblem(ode_f2, u0, tspan, param)
sol0 = solve(prob0, Tsit5())
plot(sol0, lw = 2, ylims = (0, 5))
## Wrap the above problem to a fucntion
function run_prob_2bits(;init_relax, duration,relax,signal,K,n,up,cycle)
    u0 =  Float64[i for i in rand(1:22, 17)]
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = cycle)
    param = [up,0.002,K,n,0.025,0.025,p]
    prob0 = ODEProblem(ode_f2, u0, (0,init_relax), param)
    sol0 = solve(prob0, Tsit5())
    prob = ODEProblem(ode_f2, sol0[end], tspan, param)
    sol = solve(prob, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
    return sol, ts
end

## 2 bits plot
function Bit2_plot(sol,up)
    Ylim = 3*up
    C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    B1_plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),size=(600,500),title = ["Carrying Bit with up: $up"  "1st Bit"  "2nd Bit"])
    # display(Bit2_plt)
end

# up = 1.5
# sol, ts = run_prob_2bits(;init_relax = 2000., duration=270.,relax=2000.,signal=20,K=0.031,n=1.5, up= up,cycle = 20)
# plt = Bit2_plot(sol,up)

##
anim = @animate for up = 2:0.1:3
    # up = 1.5
    sol, ts = run_prob_2bits(;init_relax = 2000., duration=270.,relax=2000.,signal=20,K=0.031,n=1.5, up= up,cycle = 20)
    plt = Bit2_plot(sol,up)
end
gif(anim, "bit2_fps15.gif", fps = 1)
##
#  ==== 2
# QS = sol[10,:]
# plot(sol, vars=(0,10))
# findpeaks(QS, m=10)#findlocalmaxima(cc.a)
#
# tsample = collect(0:sol.t[end])
# plot(tsample,sol(tsample))
# Y =collect(QS)
# filtered_Y = Y
# filtered_Y[ Y .< 1] .= 0.0
#
#
# locs = findpeaks(filtered_Y; m = 3)
# marklist = [i[1] for i in locs]
# marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
# # scatter!(sol[marklist], vars = [:g3], xlims = (1:2000.0), marker = marker_spec) # xlims = (1:t_on),
# scatter!(marklist, QS[marklist],marker = marker_spec) # xlims = (1:t_on),

# ## -------- Finding local minmax
# C_id = 10 # connector ID in the Diff equation
# opt = L_max(sol,C_id, ts[1],ts[2])
# py = plot(sol, vars=(0,C_id), plotdensity=10000);
# scatter!(py, [opt.minimizer],[-opt.minimum],label="Local Max")
##

##
# # -------- First g3 peak detection
# C_id = 10 # connector ID in the Diff equation
# opt2 = L_max(sol,C_id, ts[3],ts[4])
# opt2.minimizer, -opt2.minimum
# # --------- Attenuation averaging
# function attenuation(sol ;thred = 0.5)
#     if -opt2.minimum > thred
#         lmax_set = []
#         for i = collect(3:4:2*cycle)
#             opt_i = L_max(sol,C_id, ts[i],ts[i+1])
#             lmax = -opt_i.minimum
#             push!(lmax_set,lmax)
#         end
#         return sum(lmax_set)/length(lmax_set), std(lmax_set)
#     end
# end
#
# AT_mean, AT_var = attenuation(sol)
#
#
# function g3_peak_init(sol; ts = ts, thred = 0.5, C_id = 10)
#     opt = L_max(sol,C_id, ts[1],ts[2])
#     opt2 = L_max(sol,C_id, ts[3],ts[4])
#     lmax1 = -opt.minimum;    lmax2 = -opt2.minimum;
#     if lmax1 > thred && lmax2 < thred
#         g3_peak_i = 1
#     elseif lmax2 > thred && lmax1 < thred
#         g3_peak_i = 3
#     end
# end
#
#
# cycle_id = g3_peak_init(sol; thred = 0.5, C_id = 10)
#
#
# for ti = 1:2:length(ts)
#     @show ti, ti+1
#     opt = L_max(sol,C_id, ts[ti],ts[ti + 1])
#     lmax1 = -opt.minimum;
#     @show lmax1
# end


## === Check whether the gate switches or not ============
# function switching2(sol, f)
#     # println("End of relaxation time: ", Δ0)
#     G6 = sol(Δ0)[6];  G7 = sol(Δ0)[7]
#     G62 = sol(Δ0)[16];  G72 = sol(Δ0)[17]
#     @show G6, G7
#     @show G62, G72
#     g6ii_tot = g7ii_tot = 0
#     g6i_tot = g7i_tot = 0
#     if G6 > G7
#
#         for i = 1: 2: cycles
#             tii = Δ0 + i * Δ
#             ti = Δ0 + (i - 1) * Δ
#             # @show ti, tii
#             g6ii = sol(tii)[6];    g6i = sol(ti)[6]
#             g7ii = sol(tii)[7];    g7i = sol(ti)[7]
#             # @show g6i, g6ii
#             # @show g7i, g7ii
#
#             g6ii_tot += g6ii; g7ii_tot += g7ii
#             g6i_tot += g6i  ; g7i_tot += g7i
#         end
#     elseif G6 < G7
#         for i = 1: 2: cycles
#             tii = Δ0 + i * Δ
#             ti = Δ0 + (i - 1) * Δ
#             # @show ti, tii
#             g6ii = sol(tii)[6];    g6i = sol(ti)[6]
#             g7ii = sol(tii)[7];    g7i = sol(ti)[7]
#             # @show g6i, g6ii
#             # @show g7i, g7ii
#
#             g6ii_tot += g6ii; g7ii_tot += g7ii
#             g6i_tot += g6i  ; g7i_tot += g7i
#         end
#     end
#     cnt = length(collect(1: 2: cycles))
#     g6i_av = g6i_tot/cnt; g6ii_av = g6ii_tot/cnt
#     g7i_av = g7i_tot/cnt; g7ii_av = g7ii_tot/cnt
#     # @show g6i_av, g6ii_av
#     # @show g7i_av, g7ii_av
#     return g6i_av, g6ii_av, g7i_av, g7ii_av
# end

## Building cost function for param
# up = 1.8
# sol, ts = run_prob_2bits(;init_relax = 2000., duration=275.,relax=2000.,signal=20,K=0.011,n=1.5, up= 1.8,cycle = 20)
function Carrying_bit(sol, ts)
    ts_ID = []
    C1 = []
    for ti = 1:2:length(ts)
        push!(ts_ID, ti, ti+1)
        # connector 1
        opt = L_max(sol,10, ts[ti],ts[ti + 1])
        lmax1 = -opt.minimum;
        push!(C1,lmax1)
    end
    return ts_ID, C1
end
ts_ID, C1 = Carrying_bit(sol, ts)

function MB_T0(ts_ID, C1)
    thred = 0.1
    # Identify the T0
    # C1_B = C1 .> thred; C2_B = C2 .> thred
    T0 = []
    ts_cycle = collect(1:2:length(ts_ID))
    for i = 1:length(C1)
        if C1[i] > thred
            push!(T0,  ts_cycle[i], ts_cycle[i] + 1)
            break
        else
            continue
        end
    end
    return T0
end
# Identify the T0 for multibit counter
# T0 = MB_T0(ts_ID, C1) # if something goes wrong, this may give me 0

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
##
function Switch_cost(sol, Δ0, ts, T0; BNO = "B1")
    # G6 is HKCI, G7 is PhlF
    G6 =[];G7 =[]; G16 = []; G17 = [];
    if BNO == "B1"
        idx = [6,7]
        for sw in 0:2:14
            B1_1, B1_2 = Switch(sol,idx, Δ0, ts, T0 .+ sw)
            push!(G6, B1_1);     push!(G7, B1_2)
        end
        @show G6, G7
        return G6, G7

    elseif BNO == "B2"
        idx = [16, 17]
        for sw in 0:4:14
            B2_1, B2_2 = Switch(sol,idx, Δ0, ts, T0 .+ sw)
            push!(G16, B2_1);     push!(G17, B2_2)
        end
        @show G16, G17
        return G16, G17
    end
end

# G16, G17 = Switch_cost(sol, Δ0, ts, T0; BNO ="B2")
# G6, G7 = Switch_cost(sol, Δ0, ts, T0; BNO = "B1")
##
function osci(G6, up)
    thred = up/2
    G6b = G6 .> thred
    G6b1 = [G6b[i] for i = 1:2:length(G6b)]
    G6b2 = [G6b[i] for i = 2:2:length(G6b)]
    # @show G6b1, G6b2, G6b1 .+ G6b2
    return sum(G6b1 .+ G6b2), sum(G6b1 .* G6b2) # comp should be 4, and sum should be 0
end

Δ0 = 5000.

function cost_bit2(sol, ts, Δ0, up)
    thred = up/2
    # set T0 for multibit counter
    ts_ID, C1 = Carrying_bit(sol, ts)
    T0 = MB_T0(ts_ID, C1)
    @show T0
    if length(T0) > 0
        G6,  G7 = Switch_cost(sol, Δ0, ts, T0; BNO = "B1")
        G16, G17 = Switch_cost(sol, Δ0, ts, T0; BNO ="B2")
        # @show G6
        comp6, dot6 = osci(G6,up); comp7, dot7 = osci(G7,up);
        comp16, dot16 = osci(G16,up); comp17, dot17 = osci(G17,up);
        # @show osci(G6,up), osci(G7,up)
        # @show osci(G16,up), osci(G17,up)

        if comp6 ==4 && comp7 == 4 && dot6 == 0 && dot7 == 0 && comp16 ==2 &&comp17 ==2 && dot16 == 0 && dot17 == 0
            println("good")
            return costtot = 8
        else
            println("bad")
            return costtot = 0
        end

        b11 = G6 .> thred
        b12 = G7 .> thred
        b21 = G16 .> thred
        b22 = G17 .> thred

        # @show b11, b12
        # @show b21, b22

        b11[1] == 0 ? B1= b11 : B1= b12
        b21[1] == 0 ? B2= b21 : B2= b22

        @show B1
        # @show B2


        # cost to be converting binary to tenary rule
        costtot = 0
        for l in 1: length(B1)
            # @show B1[l]*2^0 + B2[l]*2^1 + B3[l]*2^2
            B1[l]*2^0 + B2[l]*2^1  == l ? cost = 0 : cost = 1
            costtot += cost
        end
    elseif length(T0) == 0
        costtot = 0
    end
    return costtot
end
# costtot = cost_bit2(sol, ts, up)
##

# function plot_2bit(sol, up, param)
#     plt0 = plot(sol, vars = [:g3])
#     plt1 = plot(sol, vars =[6,7])
#     plt2 = plot(sol, vars =[16,17])
#     plt = plot(plt0,plt1,plt2, layout = (3,1), ylims =(0.,3*up), title = ["Param: $param" "B1" "B2"])
#     display(plt)
# end

## Parameters Searching for  2bits counter
# ======== Sampling Parameters ===================
# Varying K, n, δ
df2 = DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[], up = Float64[])
tt = []
@time @showprogress for K = 0.011: 0.02:0.1, n = 1.5:0.1:3., δ = 250:5:350, A = 20., up = 1:0.1:3, Δ0 = 5000. #@showprogress
# @time @showprogress for K = 0.011, n = 1.5, δ = 270, A = 20., up = 1.:0.1:3. # this line is for test
    # solve DiffEq given parameters
    sol, ts = run_prob_2bits(;init_relax = Δ0, duration=δ, relax=Δ0, signal=A, K=K, n=n, up = up, cycle = 20)

    # Get cost for this parameters set
    costtot = cost_bit2(sol, ts, Δ0, up)

    println("K:$K n:$n, δ:$δ, A:$A, up:$up")
    @show costtot
    param = [K, n, δ, A, up]
    costtot == 8 ? push!(df2, param) : nothing

    # plot_2bit(sol,up, param)
    # count example
    push!(tt,1.)
    @show sum(tt)
    println("\n")

    # if sum(tt) >= 8.
    #     break
    # end
end
df2
##
CSV.write("2Bits_DB.csv", df2)
##





























## =============== 3 Bits counter case 📗 ==============
## ------ Import package and functions 🔴need to fix as 2bit does
using DifferentialEquations, ModelingToolkit
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
using Random, Base
using CSV, Optim, DataFrames
using ProgressMeter
include("functions.jl") #Latexify

## ==== Build multiple counter connectors : 3 Bits counter case 📗 =========
### Define a differential equation system
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t) g1(t) g2(t) g3(t) m2_LexA1(t) m2_IcaR(t) m2_CI1(t) m2_PsrA(t) m2_BM3RI(t) m2_HKCI(t) m2_PhlF(t) g21(t) g22(t) g23(t) g24(t) g25(t) g26(t) m3_LexA1(t) m3_IcaR(t) m3_CI1(t) m3_PsrA(t) m3_BM3RI(t) m3_HKCI(t) m3_PhlF(t)
@derivatives D'~t

eqs3 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF),
    # Connector 1 =============
    D(g1)      ~ ξ * hill(p)                   - deg(g1),
    D(g2)      ~ ξ * hill(m1_HKCI)             - deg(g2),
    D(g3)      ~ ξ * hill(g1 + g2)             - deg(g3), # g3 sserves as the input for the 2nd bit
    # Bit 2 =============
    D(m2_LexA1) ~ ξ * hill(m2_PhlF + g3)       - deg(m2_LexA1),
    D(m2_IcaR ) ~ ξ * hill(m2_LexA1 + g3)      - deg(m2_IcaR),
    D(m2_CI1  ) ~ ξ * hill(m2_LexA1 + m2_PhlF) - deg(m2_CI1),
    D(m2_PsrA ) ~ ξ * hill(m2_IcaR + m2_CI1)   - deg(m2_PsrA),
    D(m2_BM3RI) ~ ξ * hill(m2_PsrA)            - deg(m2_BM3RI),
    D(m2_HKCI ) ~ ξ * hill(m2_BM3RI + m2_PhlF) - deg(m2_HKCI),
    D(m2_PhlF ) ~ ξ * hill(m2_PsrA + m2_HKCI)  - deg(m2_PhlF),
    # Connector 2 (Two AND gates) =============
    # 1rst AND gate combines (out1,out2)
    D(g21)     ~ ξ * hill(m1_HKCI)              - deg(g21),
    D(g22)     ~ ξ * hill(m2_HKCI)             - deg(g22),
    D(g23)     ~ ξ * hill(g21 + g22)           - deg(g23), # g3 sserves as the input for the 2nd bit
    # 2nd AND gate combines (out g23,p)
    D(g24)     ~ ξ * hill(p)                   - deg(g24),
    D(g25)     ~ ξ * hill(g23)                 - deg(g25),
    D(g26)     ~ ξ * hill(g24 + g25)           - deg(g26), # g3 sserves as the input for the 2nd bit
    # Bit 3 =============
    D(m3_LexA1) ~ ξ * hill(m3_PhlF + g26)       - deg(m3_LexA1),
    D(m3_IcaR ) ~ ξ * hill(m3_LexA1 + g26)      - deg(m3_IcaR),
    D(m3_CI1  ) ~ ξ * hill(m3_LexA1 + m3_PhlF)  - deg(m3_CI1),
    D(m3_PsrA ) ~ ξ * hill(m3_IcaR + m3_CI1)    - deg(m3_PsrA),
    D(m3_BM3RI) ~ ξ * hill(m3_PsrA)             - deg(m3_BM3RI),
    D(m3_HKCI ) ~ ξ * hill(m3_BM3RI + m3_PhlF)  - deg(m3_HKCI),
    D(m3_PhlF ) ~ ξ * hill(m3_PsrA + m3_HKCI)   - deg(m3_PhlF)]
de3 = ODESystem(eqs3, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF, g1, g2, g3, m2_LexA1, m2_IcaR, m2_CI1, m2_PsrA,m2_BM3RI, m2_HKCI, m2_PhlF, g21, g22, g23, g24, g25, g26, m3_LexA1, m3_IcaR, m3_CI1, m3_PsrA, m3_BM3RI, m3_HKCI, m3_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f3 = ODEFunction(de3)
## Solve DifferentialEquations
u0 =  Float64[i for i in rand(1:22, 30)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 270, cycle = 20)
param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
prob0 = ODEProblem(ode_f3, u0, tspan, param)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)

## Plots for paper
up = 1.5
Ylim = 2*up
C_plt = plot(sol, vars = [:g3], label = "Connector 1", ylims =(0.,Ylim))
C2_plt = plot(sol, vars = [:g26], label = "Connector 2", ylims =(0.,Ylim))
B1_plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],legend = :topright, ylims =(0.,Ylim))
scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal")
B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, ylims =(0.,Ylim))
B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = :topright,  ylims =(0.,Ylim))
Bit2_plt = plot(C_plt,C2_plt, B1_plt,B2_plt,B3_plt,layout = (5,1),xtickfontsize=15,ytickfontsize=15, size=(1000,800),legend = :topright, legendfont = font("Times new roman", 8),)
##
# savefig(Bit2_plt, "3bitdynamics.png")
## optimization for parameters searching

## Wrap the above problem to a fucntion
function run_prob_3bits(;init_relax,duration,relax,signal,K,n,up,cycle)
    u0 =  Float64[i for i in rand(1:22, 30)]
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = cycle)
    param = [up,0.002,K,n,0.025,0.025,p]
    prob0 = ODEProblem(ode_f3, u0, (0,init_relax), param)
    sol0 = solve(prob0, Tsit5())
    prob = ODEProblem(ode_f3, sol0[end], tspan, param)
    sol = solve(prob, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
    return sol, ts
end

sol, ts = run_prob_3bits(;init_relax = 5000., duration=270.,relax=5000.,signal=20.,K=0.081,n=2.81, up= 1.5, cycle=20)
# sol, ts = run_prob_3bits(;duration=285.,relax=20000.,signal=10.,K=0.091,n=2.9)

## Building cost function for param

function Carrying_bit(sol, ts)
    ts_ID = []
    C1 = []; C2 = [];
    for ti = 1:2:length(ts)
        push!(ts_ID, ti, ti+1)
        # connector 1
        opt = L_max(sol,10, ts[ti],ts[ti + 1])
        lmax1 = -opt.minimum;

        # connector 2
        opt2 = L_max(sol,20, ts[ti],ts[ti + 1])
        lmax2 = -opt2.minimum;

        push!(C1,lmax1)
        push!(C2,lmax2)
    end
    return ts_ID, C1, C2
end
# ts_ID, C1, C2  = Carrying_bit(sol, ts)

# function MB_T0(ts_ID, C1, C2, up)
#     thred = 0.1
#     # Identify the T0
#     # C1_B = C1 .> thred; C2_B = C2 .> thred
#     T0 = []
#     ts_cycle = collect(1:2:length(ts_ID))
#     for i = 1:length(C1)
#         if C1[i] > thred
#             push!(T0,  ts_cycle[i], ts_cycle[i] + 1)
#             break
#         else
#             continue
#         end
#     end
#     return T0
# end
# T00 = MB_T0(ts_ID, C1, C2, up)
# Identify the T0 for multibit counter
function MB_T02(ts_ID, C1, C2, up)
    d1 = C1 .> 0.1 # return which cycle the connector 1 is above thred
    d2 = C2 .> 0.1 # return which cycle the connector 2 is above thred
    # @show d1, d2
    d1_idx = findall(x->x==1, d1)
    d2_idx = findall(x->x==1, d2)
    # @show d1_idx, d2_idx
    if isempty(d1_idx) || isempty(d2_idx) || isempty(intersect(d1_idx,d2_idx))
        println("one of the index set is empty")
        return T0 = []
    else
        cycle_id = intersect(d1_idx,d2_idx)[1]
        tid = collect(1:2:length(ts_ID))[cycle_id]
        return T0 = [tid, tid + 1 ]
    end
end
# T0 = MB_T02(ts_ID, C1, C2, up) # if something goes wrong, this may give me 0


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

function Switch_cost(sol,  Δ0, ts, T0; BNO = "B1")
    # G6 is HKCI, G7 is PhlF
    G6 =[];G7 =[]; G16 = []; G17 = [];G29 = [];G30 = [];
    if BNO == "B1"
        idx = [6,7]
        for sw in 0:2:14
            B1_1, B1_2 = Switch(sol,idx, Δ0, ts, T0 .+ sw)
            push!(G6, B1_1);     push!(G7, B1_2)
        end
        @show G6, G7
        return G6, G7

    elseif BNO == "B2"
        idx = [16, 17]
        for sw in 0:4:12
            B2_1, B2_2 = Switch(sol,idx, Δ0, ts, T0 .+ sw)
            push!(G16, B2_1);     push!(G17, B2_2)
        end
        @show G16, G17
        return G16, G17
    elseif BNO == "B3"
        idx = [29,30]
        for sw in 0:8:24
            B3_1, B3_2 = Switch(sol,idx,Δ0, ts, T0 .+ sw)
            push!(G29, B3_1);     push!(G30, B3_2)
        end
        @show G29, G30
        return G29, G30
    end
end


function osci(G6, up)
    thred = up/2
    G6b = G6 .> thred
    G6b1 = [G6b[i] for i = 1:2:length(G6b)]
    G6b2 = [G6b[i] for i = 2:2:length(G6b)]
    return sum(G6b1 .+ G6b2), sum(G6b1 .* G6b2) # comp should be 4, and sum should be 0
end


function cost_bit3(sol, ts, up)
    thred = up/2
    # set T0 for multibit counter
    ts_ID, C1, C2  = Carrying_bit(sol, ts)
    # T0 = MB_T0(ts_ID, C1, C2, up)
    T0 = MB_T02(ts_ID, C1, C2, up)
    @show T0
    if length(T0) > 0
        G6,  G7 =  Switch_cost(sol, Δ0, ts, T0; BNO = "B1")
        G16, G17 = Switch_cost(sol, Δ0, ts, T0; BNO = "B2")
        G29, G30 = Switch_cost(sol, Δ0, ts, T0; BNO = "B3")

        comp6, dot6 = osci(G6,up); comp7, dot7 = osci(G7,up);
        comp16, dot16 = osci(G16,up); comp17, dot17 = osci(G17,up);
        comp29, dot29 = osci(G29,up); comp30, dot30 = osci(G30,up);

        if comp6 ==4 && comp7 == 4 && dot6 == 0 && dot7 == 0 && comp16 ==2 &&comp17 ==2 && dot16 == 0 && dot17 == 0 && comp29 == 2 && comp30 == 2 && dot29 == 0 && dot30 == 0
            println("good")
            return costtot = 8
        else
            println("bad")
            return costtot = 0
        end

        b11 = G6 .> thred
        b12 = G7 .> thred
        b21 = G16 .> thred
        b22 = G17 .> thred
        b31 = G29 .> thred
        b32 = G30 .> thred
        # @show b11, b12
        # @show b21, b22
        # @show b31, b32
        b11[1] == 0 ? B1= b11 : B1= b12
        b21[1] == 0 ? B2= b21 : B2= b22
        b31[1] == 0 ? B3= b31 : B3= b32
        # @show B1
        # @show B2
        # @show B3

        # cost to be converting binary to tenary rule
        costtot = 0
        for l in 1: length(B1)
            # @show B1[l]*2^0 + B2[l]*2^1 + B3[l]*2^2
            B1[l]*2^0 + B2[l]*2^1 + B3[l]*2^2   == l ? cost = 0 : cost = 1
            costtot += cost
        end
    elseif length(T0) == 0
        costtot = 0
    end
    return costtot
end

# Δ0 = 20000.
# costtot = cost_bit3(sol, ts, up)
##
function plot_3bit(sol, ts, up, param)
    C_plt = plot(sol, vars = [:g3], label = "Connector 1", ylims =(0.,Ylim))
    C2_plt = plot(sol, vars = [:g26], label = "Connector 2", ylims =(0.,Ylim))
    B1_plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],legend = :topright, ylims =(0.,Ylim))
    scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal")
    B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, ylims =(0.,Ylim))
    B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = :topright,  ylims =(0.,Ylim))
    Bit3_plt = plot(C_plt,C2_plt, B1_plt,B2_plt,B3_plt,layout = (5,1),xtickfontsize=15,ytickfontsize=15, size=(1000,800),legend = :topright, legendfont = font("Times new roman", 8), title = ["C1 & Param: $param" "C2" "B1" "B2" "B3"])
    display(Bit3_plt)
end
## Parameters Searching for  3bits counter
# ======== Sampling Parameters ===================
# Varying K, n, δ
df3 = DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[], up = Float64[])
tt = []
@time @showprogress for K = 0.031: 0.02:0.1, n = 1.5:0.1:4., δ = 250:5:350, A = 20., up = 1:0.1:3 #@showprogress
    # solve DiffEq given parameters
    sol, ts = run_prob_3bits(;init_relax = 20000.,duration=δ, relax=20000., signal=A, K=K, n=n, up = up, cycle = 20)
    # Get cost for this parameters set
    costtot = cost_bit3(sol, ts, up)
    println("K:$K n:$n, δ:$δ, A:$A, up:$up\n")
    @show costtot
    param = [K, n, δ, A, up]
    costtot == 8 ? push!(df3, param) : nothing

    plot_3bit(sol,ts, up, param)
    # count example
    push!(tt,1.)
    @show sum(tt)
    # if sum(tt) >= 8.
    #     break
    # end
end
##
df3
CSV.write("3Bits_DB.csv", df3)




# ## Test df return correct set
# df = CSV.read("3bits_db.csv")
# size(df)
#
# df_correct =DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[])
# @time @showprogress for i =1:size(df)[1]
#     # i = rand(1:size(df)[1])
#     δδ = df[i,:].δ; AA =df[i,:].A;  KK = df[i,:].K; nn = df[i,:].n
#     param = [KK, nn , δδ, AA ]
#     sol, ts = run_prob_3bits(;duration=δδ, relax=20000., signal=AA, K=KK, n=nn, up =1.5)
#     ts_ID, C1, C2  = Carrying_bit(sol, ts)
#     T0 = MB_T0(ts_ID, C1, C2; thred =up/2)
#     costtot = cost_bit3(sol, ts)
#     @show costtot
#     # costtot == 8 ? push!(df_correct, param) : nothing
#     if costtot == 8
#         println(" Number $m")
#         println("δ: $δδ,  A: $AA, K: $KK,  n: $nn \n")
#         plt = plot3bits(sol,param)
#         display(plt)
#     end
# end
#
# CSV.write("3Bits_DB_corrected.csv", df_correct)

## Plots
function plot_3bit(sol, ts, up, param)
    Ylim = 3*up
    C_plt = plot(sol, vars = [:g3], label = "Connector 1", ylims =(0.,Ylim))
    C2_plt = plot(sol, vars = [:g26], label = "Connector 2", ylims =(0.,Ylim))
    B1_plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],legend = :topright, ylims =(0.,Ylim))
    scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal",ylims =(0.,Ylim))
    B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, ylims =(0.,Ylim))
    B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = :topright,  ylims =(0.,Ylim))
    Bit2_plt = plot(C_plt,C2_plt, B1_plt,B2_plt,B3_plt,layout = (5,1),xtickfontsize=15,ytickfontsize=15, size=(1000,800),legend = :topright, legendfont = font("Times new roman", 8), title = ["C1 & Param: $param" "C2" "B1" "B2" "B3"])
    display(Bit2_plt)
end

##

δδ=295.0;  AA=30.0; KK=0.091;  nn=1.8 # this is problematic
param = [δδ,  AA,  KK, nn] # outlier
@show param
sol, ts = run_prob_3bits(;duration=δδ, relax=20000., signal=AA, K=KK, n=nn)
plot3bits(sol,param)
##
costtot = cost_bit3(sol, ts)



## Check 3Bits_DB
db3 = CSV.read("/Users/chentianchi/Desktop/3bt_DB_gen/3Bits_DB.csv")
db3
for i =1:size(db3)[1]
    # i = rand(1:size(db3)[1])
    δδ = db3[i,:].δ; AA =db3[i,:].A;  KK = db3[i,:].K; nn = db3[i,:].n; upp = db3[i,:].up
    sol, ts = run_prob_3bits(;init_relax = 5000., duration=δδ,relax=5000.,signal=AA, K=KK, n=nn, up= upp, cycle=20)
    param = [KK, nn, δδ, AA, upp]
    plot_3bit(sol, ts, upp, param)
    cost_bit3(sol, ts, upp)
end
##












# ## ==== Build multiple counter connectors : 4 Bits counter case 📗 =========
# # ==== Define ODEProblem =======
# γ = 0.025
# ξ = 0.025
# hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
# degradation(x) = γ * x
# up = 1.5;
# dn = 0.002;
# K = 0.081;
# n = 2.81;
# silico4 = @ode_def_bare bit4 begin
#     # Bit 1 =================
#     dm_LexA1 = ξ * hill(m_PhlF + p) - degradation(m_LexA1)
#     dm_IcaR  = ξ * hill(m_LexA1 + p) - degradation(m_IcaR)
#     dm_CI1   = ξ * hill(m_LexA1 + m_PhlF) - degradation(m_CI1)
#     dm_PsrA  = ξ * hill(m_IcaR + m_CI1) - degradation(m_PsrA)
#     dm_BM3RI = ξ * hill(m_PsrA) - degradation(m_BM3RI)
#     dm_HKCI  = ξ * hill(m_BM3RI + m_PhlF) - degradation(m_HKCI)
#     dm_PhlF  = ξ * hill(m_PsrA + m_HKCI) - degradation(m_PhlF)
#     # Connector 1 =============
#     dg1 = ξ * hill(p) - degradation(g1)
#     dg2  = ξ * hill(m_HKCI) - degradation(g2)
#     dg3   = ξ * hill(g1 + g2) - degradation(g3) # g3 sserves as the input for the 2nd bit
#     # Bit 2 =============
#     dm2_LexA1 = ξ * hill(m2_PhlF + g3) - degradation(m2_LexA1)
#     dm2_IcaR  = ξ * hill(m2_LexA1 + g3) - degradation(m2_IcaR)
#     dm2_CI1   = ξ * hill(m2_LexA1 + m2_PhlF) - degradation(m2_CI1)
#     dm2_PsrA  = ξ * hill(m2_IcaR + m2_CI1) - degradation(m2_PsrA)
#     dm2_BM3RI = ξ * hill(m2_PsrA) - degradation(m2_BM3RI)
#     dm2_HKCI  = ξ * hill(m2_BM3RI + m2_PhlF) - degradation(m2_HKCI)
#     dm2_PhlF  = ξ * hill(m2_PsrA + m2_HKCI) - degradation(m2_PhlF)
#     # Connector 2 (Two AND gates) =============
#     # 1rst AND gate combines (out1,out2)
#     dg21 = ξ * hill(m_HKCI) - degradation(g21)
#     dg22  = ξ * hill(m2_HKCI) - degradation(g22)
#     dg23   = ξ * hill(g21 + g22) - degradation(g23) # g3 sserves as the input for the 2nd bit
#     # 2nd AND gate combines (out g23,p)
#     dg24 = ξ * hill(p) - degradation(g24)
#     dg25  = ξ * hill(g23) - degradation(g25)
#     dg26   = ξ * hill(g24 + g25) - degradation(g26) # g3 sserves as the input for the 2nd bit
#     # Bit 3 =============
#     dm3_LexA1 = ξ * hill(m3_PhlF + g26) - degradation(m3_LexA1)
#     dm3_IcaR  = ξ * hill(m3_LexA1 + g26) - degradation(m3_IcaR)
#     dm3_CI1   = ξ * hill(m3_LexA1 + m3_PhlF) - degradation(m3_CI1)
#     dm3_PsrA  = ξ * hill(m3_IcaR + m3_CI1) - degradation(m3_PsrA)
#     dm3_BM3RI = ξ * hill(m3_PsrA) - degradation(m3_BM3RI)
#     dm3_HKCI  = ξ * hill(m3_BM3RI + m3_PhlF) - degradation(m3_HKCI)
#     dm3_PhlF  = ξ * hill(m3_PsrA + m3_HKCI) - degradation(m3_PhlF)
#     # Connector 3 (Two AND gates) =============
#     # 1rst AND gate combines (out1,out2)
#     dg31 = ξ * hill(m2_HKCI) - degradation(g31)
#     dg32  = ξ * hill(m3_HKCI) - degradation(g32)
#     dg33   = ξ * hill(g31 + g32) - degradation(g33) # g3 sserves as the input for the 2nd bit
#     # 2nd AND gate combines (out g33,p)
#     dg34 = ξ * hill(p) - degradation(g34)
#     dg35  = ξ * hill(g33) - degradation(g35)
#     dg36   = ξ * hill(g34 + g35) - degradation(g36) # g3 sserves as the input for the 2nd bit
#     # Bit 4 =============
#     dm4_LexA1 = ξ * hill(m4_PhlF + g36) - degradation(m4_LexA1)
#     dm4_IcaR  = ξ * hill(m4_LexA1 + g36) - degradation(m4_IcaR)
#     dm4_CI1   = ξ * hill(m4_LexA1 + m4_PhlF) - degradation(m4_CI1)
#     dm4_PsrA  = ξ * hill(m4_IcaR + m4_CI1) - degradation(m4_PsrA)
#     dm4_BM3RI = ξ * hill(m4_PsrA) - degradation(m4_BM3RI)
#     dm4_HKCI  = ξ * hill(m4_BM3RI + m4_PhlF) - degradation(m4_HKCI)
#     dm4_PhlF  = ξ * hill(m4_PsrA + m4_HKCI) - degradation(m4_PhlF)
# end p
#
# ##
# u0 = Float64[i for i in rand(1:22, 43)]
# Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 275, cycle = 20)
# prob0 = ODEProblem(silico4, u0, tspan, p)
# sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
# Ylim = 3*up
# C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# C2_plt = plot(sol, vars = [:g26], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# C3_plt = plot(sol, vars = [:g36], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# B4_plt = plot(sol, vars = [:m4_HKCI, :m4_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
# Bit2_plt = plot(C_plt,C2_plt,C3_plt,B1_plt,B2_plt,B3_plt,B4_plt,layout = (7,1),size=(600,500))
# ##
#
# using Blink
# w = Window()
# body!(w, Bit2_plt)
