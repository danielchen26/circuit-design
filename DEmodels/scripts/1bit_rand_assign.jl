# randomly assign simulated gates for 1-bit counter.
#  Modified from silico_sim.jl 1bit section

#Author: Tianchi Chen
## ------ Import package and functions
using ModelingToolkit, OrdinaryDiffEq #DifferentialEquations
using Plots;# gr(fontfamily = "Souce Code Pro for Powerline"); #pyplot()#
using Latexify, Random, Base
using CSV, DataFrames, ProgressMeter
# include("functions.jl")# using BlackBoxOptim, LinearAlgebra
function cb_gen(ts, index, vars...)
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
        for i in eachindex(ts)
            if integrator.t == ts[i]
                integrator.p[index] = vars[i]
            end
        end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    # @show vars
    return ts, cb
end

function signal_gen(cycle, Δ0,  Δ,  δ,  A)
    # t0 = [Δ0, Δ0+ δ]; time = []; push!(time, t0[1], t0[2]);
    signal = [A, 0.]
    time = [];
    T_i = [Δ0, Δ0+ δ]
    push!(time, T_i[1], T_i[2]);
    for i in 1:cycle
        async = rand(1.:0.1:2); asyncΔ = async*Δ;
        # println("increase: ", asyncΔ)
        @. T_i += asyncΔ
        # println("time: ",T_i, diff(T_i))
        push!(time, T_i[1], T_i[2])
        push!(signal, A, 0.)
    end
    return time, signal
end

function init_control(; index = 7, Δ0 = 1000., Δ = 1000., δ = 270., cycle = 5, A = 20, p = 0.0)
    Δ0 = Δ0; Δ = Δ; δ = δ; cycle = cycle; A = A
    # async = rand(1:3); asyncΔ = async*Δ;
    # tspan = (0.0, Δ0 + cycle*asyncΔ + δ + 500.)
    time, signal = signal_gen(cycle, Δ0,  Δ,  δ, A)
    ts, cb = cb_gen([time...], index, signal...)
    p = p
    tspan = (0.0, time[end] + Δ)
    return Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p
end
# ======= find local maximum ========
function L_max(sol, var_id, ti, tf)
    f = (t) -> -sol(first(t),idxs=var_id)
    opt = optimize(f,ti,tf)
    return opt
end
## ------ Import package and functions

## ==== Build multiple counter connectors : 1Bit counter case 📗 =========
# ==== Define ODEProblem =======
# γ = 0.025
# ξ = 0.025
# up = 1.5;
# dn = 0.002;
# K = 0.081;
# n = 2.81;
# mutable struct Hill{T}
# 	dn::T
#     up::T
#     K::T
#     n::T
# end
# # @unpack dn, up, K, n = Hill(0.002, 1.5, 0.081, 2.81)
# mono_param = Hill(0.002, 1.5, 0.081, 2.81)
# hill(param::Hill, x) = param.dn + (param.up - param.dn) * param.K^(param.n) / (param.K^(param.n) + abs(x)^(param.n))


hill(dn, up, K, n, x) = dn + (up - dn) * K^(n) / (K^(n) + abs(x)^(n))
deg(x) = γ * x
#-------------------------- Define a differential equation system
@parameters t γ ξ p
@parameters LexA1[1:4] IcaR[1:4] CI1[1:4] PsrA[1:4] BM3RI[1:4] HKCI[1:4] PhlF[1:4]
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
@derivatives D'~t
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(LexA1..., m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(IcaR...,  m1_LexA1 + p)           - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(CI1...,   m1_LexA1 + m1_PhlF)         - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(PsrA...,  m1_IcaR + m1_CI1)       - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(BM3RI..., m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(HKCI...,  m1_BM3RI + m1_PhlF)     - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(PhlF...,  m1_PsrA + m1_HKCI)      - deg(m1_PhlF)]
@named de1 = ODESystem(eqs1, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF],
[ LexA1..., IcaR..., CI1..., PsrA..., BM3RI..., HKCI..., PhlF..., γ, ξ, p])
ode_f1 = ODEFunction(de1)


## ---------Run the 1bit counter problem with a single shared parameter set
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
## ----------------------------------------------------------------

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
"Remember to change the inital control parameter index, here 31 instead of 7.
31 because the 1bit counter contains 31 parameters, the control p hass the last index"
function run_prob_1bit(;init_relax, duration,relax,signal, gate_p_set::gate_param_assign)
    u0 =  rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(index = 31, Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = 3)
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

sol, ts = run_prob_1bit(;init_relax = 20000., duration=270.,relax=20000., signal=20.,K=0.081,n=2.81, up= 1.5);
up= 1.5
py = plot(sol, vars = [:m1_HKCI,:m1_PhlF],
          lw = 1.5, xlabel = "time", ylabel = "concentration",
          title = "Signal Duration: 270", ylims = (0.,3*up))
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

# cost_bit1(sol, ts, up)


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





















## ====================================== import 1bit database ==================================
using CSV, DataFrames
using Statistics, StatsPlots, DataVoyager
using Lazy:@>
using VegaLite, DataFrames, DataFramesMeta
using Distances
using LaTeXStrings
using KernelDensity, StatsBase
using Plots; pyplot()
df_1b = CSV.read("1Bit_DB_(K = 0.01: 0.05:1., n = 1.:0.2:10., δ = 10:5:500, A = 20., up = 1:1.:10).csv", DataFrame)
df_1b[2577,:]


rand_idx = rand(1:nrow(df_1b), 7)

# -----------single shared parameter casse

function gate_p_set_gen(idx, df; shared = true)
    if shared == true
        gate_p_set = gate_param_assign([0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n])
    elseif shared == "random"
        println("Will randomly select 7 index")
        rand_idx = rand(1:nrow(df), 7)
        println("Random index: ", rand_idx)
        gate_p_set = gate_param_assign([0.002, df[rand_idx[1],:].up, df[rand_idx[1],:].K, df[rand_idx[1],:].n],
                               [0.002, df[rand_idx[2],:].up, df[rand_idx[2],:].K, df[rand_idx[2],:].n],
                               [0.002, df[rand_idx[3],:].up, df[rand_idx[3],:].K, df[rand_idx[3],:].n],
                               [0.002, df[rand_idx[4],:].up, df[rand_idx[4],:].K, df[rand_idx[4],:].n],
                               [0.002, df[rand_idx[5],:].up, df[rand_idx[5],:].K, df[rand_idx[5],:].n],
                               [0.002, df[rand_idx[6],:].up, df[rand_idx[6],:].K, df[rand_idx[6],:].n],
                               [0.002, df[rand_idx[7],:].up, df[rand_idx[7],:].K, df[rand_idx[7],:].n])
    elseif shared == "gaussian" # the name can be changed, gaussion means the δ for the 7 parameter set should be close
        println("Make sure the given idx should be an array of length 7")
        gate_p_set = gate_param_assign([0.002, df[idx[1],:].up, df[idx[1],:].K, df[idx[1],:].n],
                               [0.002, df[idx[2],:].up, df[idx[2],:].K, df[idx[2],:].n],
                               [0.002, df[idx[3],:].up, df[idx[3],:].K, df[idx[3],:].n],
                               [0.002, df[idx[4],:].up, df[idx[4],:].K, df[idx[4],:].n],
                               [0.002, df[idx[5],:].up, df[idx[5],:].K, df[idx[5],:].n],
                               [0.002, df[idx[6],:].up, df[idx[6],:].K, df[idx[6],:].n],
                               [0.002, df[idx[7],:].up, df[idx[7],:].K, df[idx[7],:].n])
    end
    return gate_p_set
end


# set the shared single parameter set index
idx = 1000
idx_set = [1,2,3,4,5,6,7]
gate_p_set = gate_p_set_gen(idx_set, dff, shared="gaussian")
sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
plot(sol, vars = [m1_HKCI, m1_PhlF],label =["Q" L"\overline{Q}"])




# TO_DO
# 1. Write a gaussion filter in (n,k) space with k nearest neighour
# the (n,k) will be the 2d guassion mean, and we vary the σ to check the multual kl divergence of the 7 (n,k) δ distribution
# if the we got kl is non zero or above some threshould, we will use the mean value of the δ for the group.

# df_1b[1,:]
# first(df_1b,10)
# df_1b[df_1b.K .==[0.011,0.091],:].δ


# df_new = df_1b[df_1b.δ .== 300,:]
# idx_set = [423,430,433,440,450,451,452]
# df_new[idx_set,:]


# median(idx_set)

# df_ecoli = CSV.read("DEmodels/param_db/para_s4.csv")
# # v = Voyager(df_ecoli)



# gd = groupby(df_1b, [:n,:K])
# @transform(gd, mean_δ = mean(:δ), var_δ = std(:δ))
# using Distances
# evaluate(KLDivergence(),rand(4), rand(4))

# using Distributions
# dist = Normal(3, 1)
# cdf(dist, 1)















# dff = @> begin
#     df_1b
#     @where(:δ .∈ Ref([300.]))
#     # @by(:δ, mean_δ = mean(:δ), std_δ = std(:δ))
# end

# for i = 1:10
#     idx_set = rand(1:nrow(dff), 7)
#     dff[idx_set,:]
#     local gate_p_set = gate_p_set_gen(idx_set, dff, shared="gaussian")
#     local median_δ = dff[Int64(median(idx_set)),:].δ
#     local sol, ts = run_prob_1bit(;init_relax = 5000., duration = median_δ, relax=5000., signal=20., gate_p_set);
#     p = plot(sol, vars = [:m1_HKCI, :m1_PhlF])
#     display(p)
# end




# @where(df_1b, :K .==0.96, :n .==9.2)
# @where(df_1b, @. (:K .== 0.96, :n .== 9.2) | (:K .== 0.26,:n .== 4.0)) # not working
# @where(df_1b, @. (:K == 0.96) & (:n == 9.2) | (:K == 0.26) & (:n == 4.0)) # works


## ===================== Run from here =====================
## ===================== Run from here =====================
## ===================== Run from here =====================



"Generate a dataframe which has 7 randomly selected gates conditioned on a paticular δ"
function df_7rand_gen(df_1b ; δ = 300)
    dff = @subset(df_1b, :δ .== δ) # subset of specific δ
    idx_set = rand(1:nrow(dff), 7) # 7 random row indexes
    dff_rand = sort!(dff[idx_set,:])
    return dff, dff_rand, idx_set
end

"Generate random selected 7 gates; the subset of gates contains same n,K in the database; the 7 sample gates' index in the database"
function df_7rand_δ_dist(; style = "histogram")
    # ---generate 1 example of 7 random gate condition on default value δ = 300
    dff, dff_rand, idx_set = df_7rand_gen(df_1b)
    # selct a subset of the parameters from N-bit database that contains (n,K) in dff_rand.
    dff_δ = @> begin
        df_1b
        # @where((:K .==0.96, :n .==9.2) .| (:K .==0.26, :n .==4.0))
        @subset(@. (:K == dff_rand.K[1]) & (:n == dff_rand.n[1]) | (:K == dff_rand.K[2]) & (:n == dff_rand.n[2]) |
                (:K == dff_rand.K[3]) & (:n == dff_rand.n[3]) | (:K == dff_rand.K[4]) & (:n == dff_rand.n[4]) |
                (:K == dff_rand.K[5]) & (:n == dff_rand.n[5]) | (:K == dff_rand.K[6]) & (:n == dff_rand.n[6]) |
                (:K == dff_rand.K[7]) & (:n == dff_rand.n[7]))
        # showall()
    end
    string.(zip(dff_δ.K, dff_δ.n))
    insertcols!(dff_δ,       # DataFrame to be changed
        1,                # insert as column 1
        :point => string.(zip(dff_δ.K, dff_δ.n)),   # populate as "Day" with 1,2,3,..
    )
    if style == "boxplot"
        plt = dff_δ |>
        @vlplot(
                    width=500,
                    height=200,
                    mark={:boxplot, extent="min-max"},
                    x="point:o",
                    y={:δ, axis={title="Duration δ"}},
                ) |> display
    elseif style == "histogram"
        plt = dff_δ |>
        @vlplot(
            width=500,
            height=60,
            :bar,
            x={:δ, bin={binned=true,step=10}},
            y={"count()",title = "Count"},
            color ={:point, title = "Gates (K, n)"},
            row ={:point, title = "7 Gates δ distribution"},
            # Below config is to adjust the fontsize of everything in the grouped plot 
            # config = {
            #     axis = {titleFontSize = 18, labelFontSize = 15},
            #     legend ={titleFontSize = 25, padding =0, labelFontSize = 25},
            #     header ={labelFontSize= 15, titleFontSize = 25}
            #      }
            )  |> save("./7gates_plots/hist_delta.svg")
    end
    return dff, dff_rand, idx_set, dff_δ
end
dff, dff_rand, idx_set, dff_δ = df_7rand_δ_dist()
# pyplot()

function gate_p_set_gen(idx, df; shared = true)
    if shared == true
        gate_p_set = gate_param_assign([0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n],
                               [0.002, df[idx,:].up, df[idx,:].K, df[idx,:].n])
    elseif shared == "random"
        println("Will randomly select 7 index")
        rand_idx = rand(1:nrow(df), 7)
        println("Random index: ", rand_idx)
        gate_p_set = gate_param_assign([0.002, df[rand_idx[1],:].up, df[rand_idx[1],:].K, df[rand_idx[1],:].n],
                               [0.002, df[rand_idx[2],:].up, df[rand_idx[2],:].K, df[rand_idx[2],:].n],
                               [0.002, df[rand_idx[3],:].up, df[rand_idx[3],:].K, df[rand_idx[3],:].n],
                               [0.002, df[rand_idx[4],:].up, df[rand_idx[4],:].K, df[rand_idx[4],:].n],
                               [0.002, df[rand_idx[5],:].up, df[rand_idx[5],:].K, df[rand_idx[5],:].n],
                               [0.002, df[rand_idx[6],:].up, df[rand_idx[6],:].K, df[rand_idx[6],:].n],
                               [0.002, df[rand_idx[7],:].up, df[rand_idx[7],:].K, df[rand_idx[7],:].n])
    elseif shared == "7diff" # the name can be changed, gaussion means the δ for the 7 parameter set should be close
        println("Make sure the given idx should be an array of length 7")
        gate_p_set = gate_param_assign([0.002, df[idx[1],:].up, df[idx[1],:].K, df[idx[1],:].n],
                               [0.002, df[idx[2],:].up, df[idx[2],:].K, df[idx[2],:].n],
                               [0.002, df[idx[3],:].up, df[idx[3],:].K, df[idx[3],:].n],
                               [0.002, df[idx[4],:].up, df[idx[4],:].K, df[idx[4],:].n],
                               [0.002, df[idx[5],:].up, df[idx[5],:].K, df[idx[5],:].n],
                               [0.002, df[idx[6],:].up, df[idx[6],:].K, df[idx[6],:].n],
                               [0.002, df[idx[7],:].up, df[idx[7],:].K, df[idx[7],:].n])
    end
    return gate_p_set
end

"Test if the 7 sampled random gates generate a feasible counter"
function test()
    local gate_p_set = gate_p_set_gen(idx_set, dff, shared="7diff")
    local median_δ = dff[Int64(median(idx_set)),:].δ
    local sol, ts = run_prob_1bit(;init_relax = 5000., duration = median_δ, relax=5000., signal=20., gate_p_set);
    p = plot(sol, vars = [m1_HKCI, m1_PhlF],
            label =["Q: output" L"$\bar{Q}$: variable to feed back"],
            xlabel = " Time steps", ylabel = "Concentration",
            xtickfontsize=15, ytickfontsize=15,
            legendfontsize= 10, guidefont=15,
            ylim = maximum(dff.up)*1.2,
            dpi = 300)
    display(p)
    return p
end

# ==== one run test =====
dff, dff_rand, idx_set, dff_δ = df_7rand_δ_dist()
p = test()
# ===== if the above two line result looks good. Then save it to run the line in below
p |> save("./7gates_plots/1bit_exmple_Lfont.png") 


# In  df_7rand_δ_dist function, I added config section to adjust the font size of everything. This is already integrated
#  plt = dff_δ |>
        # @vlplot(
        #     width=500,
        #     height=60,
        #     :bar,
        #     x={:δ, bin={binned=true,step=10}},
        #     y={"count()",title = "Count"},
        #     color ={:point, title = "Gates (K, n)"},
        #     row ={:point, title = "7 Gates δ distribution", titleFontSize = 30},
        #     # title = "test",
        #     config = {
        #         axis = {titleFontSize = 18, labelFontSize = 15},
        #         legend ={titleFontSize = 25, padding =0, labelFontSize = 25},
        #         header ={labelFontSize= 15, titleFontSize = 25}
        #          }
        #     )






## ===================== Run end from here =====================
## ===================== Run end from here =====================
## ===================== Run end from here =====================




# calculate pairwise KL distances between each column of δ distribution of 7 gates.
"Ouput the pairwise KL distances and overlap of 7 quantile numbers associated with 7 gates_δ_distribution"
function KL_δ(dff_δ; vis = true )
    gd = groupby(dff_δ,:point)
    keys(gd)
    # for each group I will sample the distribution and resample the data
    # KDE for ith point δ distribution
    U_set  = [kde(gd[i].δ) for i in 1:length(gd)]
    # resample from the estimated distribution to get 1000 data points.
    println("Ssample size from distribution: ", [size(i.density) for i in U_set])
    U_set_resample = [sample(U.x, weights(U.density), size(U.density)) for U in U_set ] # with size(U.density) number of points
    # make 7 arrary to matrix
    X = reshape(U_set_resample, :,7)
    # calculate pairwise KL divergence.
    R = pairwise(KLDivergence(), hcat(X...), dims=2)
    # calculate quantile for each δ distribution.
    qt = [quantile(sample(data.x, weights(data.density), size(data.density)),  weights(data.density),  [0.25, 0.75]) for data in U_set]
    qt1_max = maximum(first.(qt))
    qt2_min = minimum([i[2] for i in qt])
    if qt1_max < qt2_min
        println("We should obtain a feasible counter")
    end
    if vis ==true
        display(heatmap(R, color = :viridis))
        return R, qt, qt1_max, qt2_min
    end
    return R, qt, qt1_max, qt2_min
end
# one instance
R, qt, qt1_max, qt2_min = KL_δ(dff_δ)
qt
qt1_max, qt2_min
# calculate the quantile of a distribution

# qt_max should < qt_min
# qt1_max = maximum(first.(qt))
# qt2_min = minimum([i[2] for i in qt])
# if qt1_max < qt2_min
#     println("We should obtain a feasible counter")
# end

# function overlap(a,b)
#     @show c1 = a[1], b[1]
#     @show c2 = a[2], b[2]
#     return [maximum(c1),minimum(c2)]
# end

# qi = overlap(qt[1],qt[2])
# for i in 1: length(qt)
#     qi = overlap(qt[i], qt[i+1])
#     @show qi
# end





# -----give a test δ value to begin ----
"This function generate a list of dataframes which correspond to 7 randomly selected point in n,K space conditioned on a fixed δ"
function df_rand_set_gen(df_1b, ;δ =300)
    dff = @where(df_1b, :δ .== δ)
    idx_set = rand(1:nrow(dff), 7)
    dff_rand = dff[idx_set,:]

    rand_select_δ = []
    # @time [push!(rand_select_δ, @where(df_1b, :K .==i.K, :n .==i.n).δ) for i in eachrow(dff_rand) ]
    @time for i  in eachrow(dff_rand)
        df_i = @where(df_1b, :K .==i.K, :n .==i.n)
        push!(rand_select_δ, df_i.δ)
    end

    @time point_names = [string(collect(row)) for row in eachrow(dff_rand[:,[:n,:K]])]
    rand_select_δ

    df_rand_set = []
    for i = 1:7
        df_i = DataFrame(value = rand_select_δ[i], point = point_names[i])
        push!(df_rand_set, df_i)
    end
    return df_rand_set
end
df_rand_set = df_rand_set_gen(df_1b, ;δ =300)


# boxplot
"Visualization of distribution of δ values in each randomly selected gate. Style option: boxplot, histogram, violin, density"
function gates_δ_distribution(df_rand_set; style = "boxplot")
    if style == "boxplot"
        vcat(df_rand_set...) |>
        @vlplot(
            width=500,
            height=200,
            mark={:boxplot, extent="min-max"},
            x="point:o",
            y={:value, axis={title="Duration δ"}},
        )
    elseif style == "histogram"
        vcat(df_rand_set...) |>
        @vlplot(
            width=500,
            height=60,
            :bar,
            x={:value, bin={binned=true,step=10}},
            y={"count()",title = "Count"},
            color =:point,
            row =:point
        )
    elseif style == "violin"
        vcat(df_rand_set...) |>
        @vlplot(
            mark={:area, orient="horizontal"},
            transform=[
                {density="value", groupby=["point"],
                as=["value", "density"]}
            ],
            y="value:q",
            x= {"density:q", stack="center", impute=nothing, title=nothing,
                axis={labels=false, values=[0], grid=false, ticks=true}},
            column={"point:n", header={titleOrient="bottom", labelOrient="bottom",
                    labelPadding=0}},
            color = "point:n",
            width=70,
            spacing=0,
            config={view={stroke=nothing}}
        )
    elseif style == "density"
        vcat(df_rand_set...) |>
        @vlplot(
            width=500,
            height=50,
            :area,
            transform=[
                {density="value",bandwidth = 1, groupby=["point"],counts=true}
            ],
            x={"value:q", title="Duration δ"},
            y= {"density:q",stack=true},
            color={"point:n",scale={scheme=:category20}},
            opacity={value=0.8},
            row = :point
        )
    end
end




dff, dff_rand, idx_set, dff_δ = df_7rand_δ_dist()
df_rand_set = df_rand_set_gen(df_1b, ;δ =300)
# example of 4 different plots for a given δ randomly generated 7 gates.
plt_box = gates_δ_distribution(df_rand_set, style ="boxplot")
plt_hist = gates_δ_distribution(df_rand_set, style ="histogram")
plt_violin = gates_δ_distribution(df_rand_set, style ="violin")
plt_density = gates_δ_distribution(df_rand_set, style ="density")






function nK_Knn(point; direction = "++", radius = 3, knn = 3)
    @show point

end

nK_Knn([3.6,0.41])

