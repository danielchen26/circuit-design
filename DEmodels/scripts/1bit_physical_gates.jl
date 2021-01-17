## Load package
using OrdinaryDiffEq, ModelingToolkit, CSV, Statistics
using DataFramesMeta, DataFrames, LaTeXStrings, Plots;pyplot()
using StatsBase
mutable struct gate_param_assign{T}
	g1::T
    g2::T
    g3::T
    g4::T
    g5::T
    g6::T
    g7::T
end
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
function init_control(; index = 31, Δ0 = 1000., Δ = 1000., δ = 270., cycle = 5, A = 20, p = 0.0)
    Δ0 = Δ0; Δ = Δ; δ = δ; cycle = cycle; A = A
    # async = rand(1:3); asyncΔ = async*Δ;
    # tspan = (0.0, Δ0 + cycle*asyncΔ + δ + 500.)
    time, signal = signal_gen(cycle, Δ0,  Δ,  δ, A)
    ts, cb = cb_gen([time...], index, signal...)
    p = p
    tspan = (0.0, time[end] + Δ)
    return Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p
end
# load database for 1bit counter
pwd()
db = CSV.read("1Bit_DB_(K = 0.01: 0.05:1., n = 1.:0.2:10., δ = 10:5:500, A = 20., up = 1:1.:10).csv",DataFrame)
# load para_s4 table
df = CSV.read("./DEmodels/param_db/para_s4.csv",DataFrame)
rename!(df, [:repressor, :RBS, :dn, :up, :K, :n], makeunique=true)
## Define model
hill(dn, up, K, n, x) = dn + (up - dn) * K^(n) / (K^(n) + abs(x)^(n))
deg(x) = γ * x

@parameters t γ ξ p
@parameters LexA1[1:4] IcaR[1:4] CI1[1:4] PsrA[1:4] BM3RI[1:4] HKCI[1:4] PhlF[1:4]
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
@derivatives D'~t
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(LexA1..., m1_PhlF + p)            - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(IcaR...,  m1_LexA1 + p)           - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(CI1...,   m1_LexA1 + m1_PhlF)     - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(PsrA...,  m1_IcaR + m1_CI1)       - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(BM3RI..., m1_PsrA)                - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(HKCI...,  m1_BM3RI + m1_PhlF)     - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(PhlF...,  m1_PsrA + m1_HKCI)      - deg(m1_PhlF)]
de1 = ODESystem(eqs1, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF],
[ LexA1..., IcaR..., CI1..., PsrA..., BM3RI..., HKCI..., PhlF..., γ, ξ, p])
ode_f1 = ODEFunction(de1)

## Generate gate parameters from library (E.coli & Yeast)
function unique_gate_sample(df)
	grouped = groupby(df,:repressor)
	g7 = grouped[sample(1:grouped.ngroups,7,replace = false)]
	vcat([dfi[sample(1:nrow(dfi), 1, replace=false),:] for dfi in g7]...)
end
function gate_p_set_gen(idx, df; shared = true, rand_num = 7)
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
        rand_idx = rand(1:nrow(df), rand_num)
        println("Random index: ", rand_idx)
		if rand_num < 7
			@show rand_num
			p_set = []
			for i in 1:rand_num
				push!(p_set,[df[rand_idx[i],:].dn, df[rand_idx[i],:].up, df[rand_idx[i],:].K, df[rand_idx[i],:].n])
			end
			@show p_set
			@show fix_id = rand_idx[1]# i.e., once 3 random gates are selected, the rest 4 fours all defaultly use the first gates parameter.
			for i in 1:7 - rand_num
				push!(p_set,[df[fix_id,:].dn, df[fix_id,:].up, df[fix_id,:].K, df[fix_id,:].n])
			end
			@show p_set
			gate_p_set = gate_param_assign(p_set...)
			dump(gate_p_set)
			return gate_p_set, vcat(rand_idx, fix_id*ones(7-rand_num))
		elseif rand_num == 7
			# # force 7 gates all have different parameter set
			# while length(unique(rand_idx)) < 7
			# 	rand_idx = rand(1:nrow(df), 7)
			# end
			# println("The final random index: ", rand_idx)
			println("random index is not relevant for 7 unique gates case")
			df_unique = unique_gate_sample(df)
			@show df_unique
        	gate_p_set = gate_param_assign([df_unique[1,:].dn, df_unique[1,:].up, df_unique[1,:].K, df_unique[1,:].n],
                               [df_unique[2,:].dn, df_unique[2,:].up, df_unique[2,:].K, df_unique[2,:].n],
                               [df_unique[3,:].dn, df_unique[3,:].up, df_unique[3,:].K, df_unique[3,:].n],
                               [df_unique[4,:].dn, df_unique[4,:].up, df_unique[4,:].K, df_unique[4,:].n],
                               [df_unique[5,:].dn, df_unique[5,:].up, df_unique[5,:].K, df_unique[5,:].n],
                               [df_unique[6,:].dn, df_unique[6,:].up, df_unique[6,:].K, df_unique[6,:].n],
                               [df_unique[7,:].dn, df_unique[7,:].up, df_unique[7,:].K, df_unique[7,:].n])
							   return gate_p_set, df_unique
		end

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



# @show unique_gate_sample(df)


@show gate_p_set,df_unique= gate_p_set_gen(rand(), df; shared = "random")
# dump(gate_p_set)


## Run one example
"Remember to change the inital control parameter index, here 31 instead of 7.
31 because the 1bit counter contains 31 parameters, the control p hass the last index"
function run_prob_1bit(;init_relax, duration,relax,signal, gate_p_set::gate_param_assign,cycle) # cycle default =10
    u0 =  rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(index = 31, Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = cycle)
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


# ## set the shared single parameter set index 🔺
# idx = 1000
# idx_set = [1,2,3,4,5,6,7]
# gate_p_set = gate_p_set_gen(idx_set, dff, shared="gaussian")
# sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
# plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"])



## The criterion for picking a feasible counter
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
    for sw in 0:2:14 # 💚changed 8 from 14
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
function cost_bit1(sol, ts) # modified to customize up values
	g67min =  maximum([minimum(sol[6,Int64(round(length(sol)/2)):end]),minimum(sol[7,Int64(round(length(sol)/2)):end])]);
	g67max = minimum([maximum(sol[6,Int64(round(length(sol)/2)):end]),maximum(sol[7,Int64(round(length(sol)/2)):end])]);
	up = (g67min + g67max)/2
	@show up
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
##


## test random selected 3 gate from para_s4 dataset to see if counter works.
# ---------------- #
# - rand_num = 3 - #
# ---------------- #
gate_p_set, rand_idx = gate_p_set_gen(rand(), df; shared = "random", rand_num = 3)
# sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
# dump(gate_p_set)
all_7_up_set = [getfield(gate_p_set,name)[2] for name in fieldnames(gate_param_assign)];up = mean(all_7_up_set);

for dt in 200:5:500
	sol, ts = run_prob_1bit(;init_relax = 5000., duration= dt, relax=5000., signal=20., gate_p_set);
	plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$dt")
	ylims!((0.0,6))
	display(plt)
end

@show convert.(Int64,rand_idx)
@show df[convert.(Int64,rand_idx),:]



## Verify the individual parameters
# ---------------- #
# - rand_num = 3 - #
# ---------------- #
# 1. example 1
# convert.(Int64, rand_idx) = [18, 9, 17, 18, 18, 18, 18]    340 < δ < 430
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 2   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 3   │ SrpR      │ S1     │ 0.003   │ 1.3     │ 0.01    │ 2.9     │
# │ 4   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 5   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 6   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 7   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# 2. example 2
# convert.(Int64, rand_idx) = [12, 16, 5, 12, 12, 12, 12]   380 < δ < 470
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 2   │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │
# │ 3   │ BM3R1     │ B2     │ 0.005   │ 0.5     │ 0.15    │ 2.9     │
# │ 4   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 5   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 6   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 7   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │

# specific_idx = [18, 9, 17, 18, 18, 18, 18] # [12, 16, 5, 12, 12, 12, 12]
specific_idx = [18, 9, 17, 18, 18, 18, 18]
# specific_idx = [12, 16, 5, 12, 12, 12, 12]
p_set_specific = []
[push!(p_set_specific,[df[specific_idx[i],:].dn, df[specific_idx[i],:].up, df[specific_idx[i],:].K, df[specific_idx[i],:].n]) for i in 1:7]
gate_p_set_specific = gate_param_assign(p_set_specific...)
δ = 400
sol, ts = run_prob_1bit(;init_relax = 5000., duration= δ, relax=5000., signal=20., gate_p_set = gate_p_set_specific);
plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$δ",dpi = 300)
display(plt)
# check cost function
costtot = cost_bit1(sol, ts, up)

# save figure
# case 1
# savefig(plt, "./DEmodels/scripts/Paper_plots/physical_gate_cases/3diff_340 < δ < 430_2.png")





## test random selected 4 gate from para_s4 dataset to see if counter works.
# ---------------- #
# - rand_num = 4 - #
# ---------------- #
gate_p_set, rand_idx = gate_p_set_gen(rand(), df; shared = "random", rand_num = 4)
# sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
# dump(gate_p_set)
all_7_up_set = [getfield(gate_p_set,name)[2] for name in fieldnames(gate_param_assign)];up = mean(all_7_up_set);

for dt in 200:5:500
	sol, ts = run_prob_1bit(;init_relax = 5000., duration= dt, relax=5000., signal=20., gate_p_set);
	plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$dt")
	ylims!((0.0,6))
	display(plt)
end

@show convert.(Int64,rand_idx)
@show df[convert.(Int64,rand_idx),:]



## Working example
# ---------------- #
# - rand_num = 4 - #
# ---------------- #
# convert.(Int64, rand_idx) = [13, 16, 1, 3, 13, 13, 13]   375 < δ < 400
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 2   │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │
# │ 3   │ AmeR      │ F1     │ 0.2     │ 3.8     │ 0.09    │ 1.4     │
# │ 4   │ BetI      │ E1     │ 0.07    │ 3.8     │ 0.41    │ 2.4     │
# │ 5   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 6   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 7   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │

# convert.(Int64, rand_idx) = [15, 15, 12, 6, 20, 15, 15]
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 2   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 3   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 4   │ BM3R1     │ B3     │ 0.01    │ 0.8     │ 0.26    │ 3.4     │
# │ 5   │ SrpR      │ S4     │ 0.007   │ 2.1     │ 0.1     │ 2.8     │
# │ 6   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 7   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │

specific_idx = [13, 16, 1, 3, 13, 13, 13]
p_set_specific = []
[push!(p_set_specific,[df[specific_idx[i],:].dn, df[specific_idx[i],:].up, df[specific_idx[i],:].K, df[specific_idx[i],:].n]) for i in 1:7]
gate_p_set_specific = gate_param_assign(p_set_specific...)
δ = 390
sol, ts = run_prob_1bit(;init_relax = 5000., duration= δ, relax=5000., signal=20., gate_p_set = gate_p_set_specific);
plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$δ", dpi = 300)
display(plt)
# check cost function
costtot = cost_bit1(sol, ts, up)

# save figure
# case 1
savefig(plt, "./DEmodels/scripts/Paper_plots/physical_gate_cases/4diff_ δ=390_2.png")


## test random selected 5 gate from para_s4 dataset to see if counter works.
# ---------------- #
# - rand_num = 5 - #
# ---------------- #
gate_p_set, rand_idx = gate_p_set_gen(rand(), df; shared = "random", rand_num = 5)
# sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
# dump(gate_p_set)
all_7_up_set = [getfield(gate_p_set,name)[2] for name in fieldnames(gate_param_assign)];up = mean(all_7_up_set);

for dt in 200:5:500
	sol, ts = run_prob_1bit(;init_relax = 5000., duration= dt, relax=5000., signal=20., gate_p_set);
	plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$dt")
	# ylims!((0.0,6))
	display(plt)
end

@show convert.(Int64,rand_idx)
@show df[convert.(Int64,rand_idx),:]

## 🔺special case for rand_num = 5 to watch out, seems to have a really wide δ range, but the steady state upper bound seems to be quite small
# # convert.(Int64, rand_idx) = [9, 9, 8, 8, 13, 9, 9]
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 2   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 3   │ IcaRA     │ I1     │ 0.08    │ 2.2     │ 0.1     │ 1.4     │
# │ 4   │ IcaRA     │ I1     │ 0.08    │ 2.2     │ 0.1     │ 1.4     │
# │ 5   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 6   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 7   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │

## Working example
# ---------------- #
# - rand_num = 5 - #
# ---------------- #
# convert.(Int64, rand_idx) = [13, 16, 15, 19, 12, 13, 13]      290 <  δ < 445
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 2   │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │
# │ 3   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 4   │ SrpR      │ S3     │ 0.004   │ 2.1     │ 0.06    │ 2.8     │
# │ 5   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 6   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 7   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
specific_idx = [13, 16, 15, 19, 12, 13, 13]
p_set_specific = []
[push!(p_set_specific,[df[specific_idx[i],:].dn, df[specific_idx[i],:].up, df[specific_idx[i],:].K, df[specific_idx[i],:].n]) for i in 1:7]
gate_p_set_specific = gate_param_assign(p_set_specific...)
δ = 400
sol, ts = run_prob_1bit(;init_relax = 5000., duration= δ, relax=5000., signal=20., gate_p_set = gate_p_set_specific);
plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$δ")
display(plt)
# check cost function
costtot = cost_bit1(sol, ts, up)


# save figure
# case 1
savefig(plt, "./DEmodels/scripts/Paper_plots/physical_gate_cases/5diff_290<δ<445_2.png")


## test random selected 6 gate from para_s4 dataset to see if counter works.
# ---------------- #
# - rand_num = 6 - #
# ---------------- #
gate_p_set, rand_idx = gate_p_set_gen(rand(), df; shared = "random", rand_num = 6)
# sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
# dump(gate_p_set)
all_7_up_set = [getfield(gate_p_set,name)[2] for name in fieldnames(gate_param_assign)];up = mean(all_7_up_set);

for dt in 200:5:500
	sol, ts = run_prob_1bit(;init_relax = 5000., duration= dt, relax=5000., signal=20., gate_p_set);
	plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$dt")
	# ylims!((0.0,6))
	display(plt)
end

@show convert.(Int64,rand_idx)
@show df[convert.(Int64,rand_idx),:]




## Working example
# ---------------- #
# - rand_num = 6 - #
# ---------------- #
# convert.(Int64, rand_idx) = [13, 7, 4, 15, 5, 2, 13]      365 < δ < 530
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 2   │ HlyIIR    │ H1     │ 0.07    │ 2.5     │ 0.19    │ 2.6     │
# │ 3   │ BM3R1     │ B1     │ 0.004   │ 0.5     │ 0.04    │ 3.4     │
# │ 4   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 5   │ BM3R1     │ B2     │ 0.005   │ 0.5     │ 0.15    │ 2.9     │
# │ 6   │ AmtR      │ A1     │ 0.06    │ 3.8     │ 0.07    │ 1.6     │
# │ 7   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │

# convert.(Int64, rand_idx) = [12, 5, 9, 6, 20, 7, 12]   415 < δ < 440
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 2   │ BM3R1     │ B2     │ 0.005   │ 0.5     │ 0.15    │ 2.9     │
# │ 3   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 4   │ BM3R1     │ B3     │ 0.01    │ 0.8     │ 0.26    │ 3.4     │
# │ 5   │ SrpR      │ S4     │ 0.007   │ 2.1     │ 0.1     │ 2.8     │
# │ 6   │ HlyIIR    │ H1     │ 0.07    │ 2.5     │ 0.19    │ 2.6     │
# │ 7   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │

# very interesting case 1.1 ! 💚💚💚💚💚  δ seems can be extremely large , δ > 4000 💚💚💚💚💚
# convert.(Int64, rand_idx) = [5, 3, 17, 12, 12, 9, 8]
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ BM3R1     │ B2     │ 0.005   │ 0.5     │ 0.15    │ 2.9     │
# │ 2   │ BetI      │ E1     │ 0.07    │ 3.8     │ 0.41    │ 2.4     │
# │ 3   │ SrpR      │ S1     │ 0.003   │ 1.3     │ 0.01    │ 2.9     │
# │ 4   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 5   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 6   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 7   │ IcaRA     │ I1     │ 0.08    │ 2.2     │ 0.1     │ 1.4     │
# very interesting case 1.2 ! 💚💚💚💚💚  δ seems can be extremely large , δ > 4000 💚💚💚💚💚
# convert.(Int64, rand_idx) = [18, 17, 1, 16, 19, 10, 10]
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 2   │ SrpR      │ S1     │ 0.003   │ 1.3     │ 0.01    │ 2.9     │
# │ 3   │ AmeR      │ F1     │ 0.2     │ 3.8     │ 0.09    │ 1.4     │
# │ 4   │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │
# │ 5   │ SrpR      │ S3     │ 0.004   │ 2.1     │ 0.06    │ 2.8     │
# │ 6   │ LmrA      │ N1     │ 0.2     │ 2.2     │ 0.18    │ 2.1     │
# │ 7   │ LmrA      │ N1     │ 0.2     │ 2.2     │ 0.18    │ 2.1     │

# very interesting case 2. ! 💚💚💚💚💚 δ has two feasible ranges
# convert.(Int64, rand_idx) = [7, 12, 7, 13, 18, 17, 6]   310  < δ  < 355    485 < δ < 545
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ HlyIIR    │ H1     │ 0.07    │ 2.5     │ 0.19    │ 2.6     │
# │ 2   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 3   │ HlyIIR    │ H1     │ 0.07    │ 2.5     │ 0.19    │ 2.6     │
# │ 4   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 5   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 6   │ SrpR      │ S1     │ 0.003   │ 1.3     │ 0.01    │ 2.9     │
# │ 7   │ BM3R1     │ B3     │ 0.01    │ 0.8     │ 0.26    │ 3.4     │

# specific_idx = [13, 7, 4, 15, 5, 2, 13]
specific_idx = [5, 3, 17, 12, 12, 9, 8]
p_set_specific = []
[push!(p_set_specific,[df[specific_idx[i],:].dn, df[specific_idx[i],:].up, df[specific_idx[i],:].K, df[specific_idx[i],:].n]) for i in 1:7]
gate_p_set_specific = gate_param_assign(p_set_specific...)
δ = 10000#530
sol, ts = run_prob_1bit(;init_relax = 5000., duration= δ, relax=5000., signal=20., gate_p_set = gate_p_set_specific);
plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$δ",dpi = 300)
display(plt)
# check cost function
costtot = cost_bit1(sol, ts, up)


# save figure
# case 1
# savefig(plt, "./DEmodels/scripts/Paper_plots/physical_gate_cases/6diff_δ ~4000_1.png")




## test random selected 7 gate from para_s4 dataset to see if counter works.
# ---------------- #
# - rand_num = 7 - #
# ---------------- #
gate_p_set,df_unique = gate_p_set_gen(rand(), df; shared = "random", rand_num = 7)
all_7_up_set = [getfield(gate_p_set,name)[2] for name in fieldnames(gate_param_assign)];up = mean(all_7_up_set);
all_7_up_set

for dt in 200:5:700
	sol, ts = run_prob_1bit(;init_relax = 5000., duration= dt, relax=5000., signal=20., gate_p_set,cycle = 5);
	# costtot = cost_bit1(sol, ts, up)
	# @show costtot, dt
	# costtot == 0 ? push!(rand_idx_δ_set, [rand_idx,dt]) && push!(gates_p_set,gate_p_set) : nothing
	plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$dt")
	display(plt)
end

p_unique = [];[push!(p_unique,[i.dn, i.up, i.K, i.n]) for i in eachrow(df_unique[:, [:dn,:up,:K, :n]])]
@show p_unique


## Working example
# ---------------- #
# - rand_num = 7 - #
# ---------------- #
# very interesting case 1.1 ! 💚💚💚💚💚  δ seems can be extremely large , 💚💚💚💚💚 δ ~ 4000
# convert.(Int64, rand_idx) = [20, 7, 19, 8, 10, 6, 5]  🍏additional rand_idx = [13, 17, 3, 18, 15, 1, 12]
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ SrpR      │ S4     │ 0.007   │ 2.1     │ 0.1     │ 2.8     │
# │ 2   │ HlyIIR    │ H1     │ 0.07    │ 2.5     │ 0.19    │ 2.6     │
# │ 3   │ SrpR      │ S3     │ 0.004   │ 2.1     │ 0.06    │ 2.8     │
# │ 4   │ IcaRA     │ I1     │ 0.08    │ 2.2     │ 0.1     │ 1.4     │
# │ 5   │ LmrA      │ N1     │ 0.2     │ 2.2     │ 0.18    │ 2.1     │
# │ 6   │ BM3R1     │ B3     │ 0.01    │ 0.8     │ 0.26    │ 3.4     │
# │ 7   │ BM3R1     │ B2     │ 0.005   │ 0.5     │ 0.15    │ 2.9     │

# very interesting case 1.2 ! 💚💚💚💚💚  δ seems can be extremely large , 💚💚💚💚💚 δ ~ 10000
# convert.(Int64, rand_idx) = [13, 17, 3, 18, 15, 1, 12]
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ PhIF      │ P3     │ 0.02    │ 6.8     │ 0.23    │ 4.2     │
# │ 2   │ SrpR      │ S1     │ 0.003   │ 1.3     │ 0.01    │ 2.9     │
# │ 3   │ BetI      │ E1     │ 0.07    │ 3.8     │ 0.41    │ 2.4     │
# │ 4   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 5   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 6   │ AmeR      │ F1     │ 0.2     │ 3.8     │ 0.09    │ 1.4     │
# │ 7   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │

# very interesting case 2.1 ! 💚💚💚💚💚  δ can have > 1 feasible range
# convert.(Int64, rand_idx) = [2, 6, 15, 20, 5, 11, 12]   355 < δ < 395     625 < δ < 650
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ AmtR      │ A1     │ 0.06    │ 3.8     │ 0.07    │ 1.6     │
# │ 2   │ BM3R1     │ B3     │ 0.01    │ 0.8     │ 0.26    │ 3.4     │
# │ 3   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 4   │ SrpR      │ S4     │ 0.007   │ 2.1     │ 0.1     │ 2.8     │
# │ 5   │ BM3R1     │ B2     │ 0.005   │ 0.5     │ 0.15    │ 2.9     │
# │ 6   │ PhIF      │ P1     │ 0.01    │ 3.9     │ 0.03    │ 4.0     │
# │ 7   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# very interesting case 2.2 ! 💚💚💚💚💚  δ can have > 1 feasible range
# convert.(Int64, rand_idx) = [6, 18, 19, 12, 15, 11, 16]     225 < δ < 395     470 < δ < 515    680 < δ < 720
# df[convert.(Int64, rand_idx), :] = 7×6 DataFrame
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ BM3R1     │ B3     │ 0.01    │ 0.8     │ 0.26    │ 3.4     │
# │ 2   │ SrpR      │ S2     │ 0.003   │ 2.1     │ 0.04    │ 2.6     │
# │ 3   │ SrpR      │ S3     │ 0.004   │ 2.1     │ 0.06    │ 2.8     │
# │ 4   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 5   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 6   │ PhIF      │ P1     │ 0.01    │ 3.9     │ 0.03    │ 4.0     │
# │ 7   │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │

# 7gates example with new sampling function
# example 1.
# │ Row │ repressor │ RBS    │ dn      │ up      │ K       │ n       │
# │     │ String    │ String │ Float64 │ Float64 │ Float64 │ Float64 │
# ├─────┼───────────┼────────┼─────────┼─────────┼─────────┼─────────┤
# │ 1   │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 2   │ PhIF      │ P2     │ 0.02    │ 4.1     │ 0.13    │ 3.9     │
# │ 3   │ LmrA      │ N1     │ 0.2     │ 2.2     │ 0.18    │ 2.1     │
# │ 4   │ PsrA      │ R1     │ 0.2     │ 5.9     │ 0.19    │ 1.8     │
# │ 5   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 6   │ AmtR      │ A1     │ 0.06    │ 3.8     │ 0.07    │ 1.6     │
# │ 7   │ HlyIIR    │ H1     │ 0.07    │ 2.5     │ 0.19    │ 2.6     │
# p_unique = Any[[0.01, 2.4, 0.05, 2.7], [0.02, 4.1, 0.13, 3.9], [0.2, 2.2, 0.18, 2.1], [0.2, 5.9, 0.19, 1.8], [0.07, 4.3, 0.05, 1.7], [0.06, 3.8, 0.07, 1.6], [0.07, 2.5, 0.19, 2.6]]
# p_unique = Any[[0.01, 3.9, 0.03, 4.0], [0.07, 4.3, 0.05, 1.7], [0.06, 3.8, 0.07, 1.6], [0.07, 3.8, 0.41, 2.4], [0.07, 2.5, 0.19, 2.6], [0.01, 2.4, 0.05, 2.7], [0.004, 2.1, 0.06, 2.8]]
# p_unique = Any[[0.003, 1.3, 0.01, 2.9], [0.2, 3.8, 0.09, 1.4], [0.07, 4.3, 0.05, 1.7], [0.07, 3.8, 0.41, 2.4], [0.2, 5.9, 0.19, 1.8], [0.02, 4.1, 0.13, 3.9], [0.2, 2.2, 0.18, 2.1]]
# p_unique = Any[[0.07, 4.3, 0.05, 1.7], [0.2, 5.9, 0.19, 1.8], [0.03, 2.8, 0.21, 2.4], [0.2, 2.2, 0.18, 2.1], [0.07, 2.5, 0.19, 2.6], [0.07, 3.8, 0.41, 2.4], [0.06, 3.8, 0.07, 1.6]]
# below is cello parameter
p_unique = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]

function check_δ_range(p_unique, δ_set)
	gate_p_unique = gate_param_assign(p_unique...)
	for δ in δ_set
		sol, ts = run_prob_1bit(;init_relax = 5000., duration= δ, relax=5000., signal=20., gate_p_set = gate_p_unique,cycle =10);
		# g67min =  maximum([minimum(sol[6,Int64(round(length(sol)/2)):end]),minimum(sol[7,Int64(round(length(sol)/2)):end])]);
		# g67max = minimum([maximum(sol[6,Int64(round(length(sol)/2)):end]),maximum(sol[7,Int64(round(length(sol)/2)):end])]);
		# up = (g67min + g67max)/2
		# @show up
		costtot = cost_bit1(sol, ts)
		plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"],legendfontsize = 3, title = L"\delta=\ %$δ")
		display(plt)
		if length(δ_set) ==1
			return plt
		else
			nothing
		end
	end

end

Hill(dn, up, K, n, x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
function plot_activation(xrange,gate)
	x = collect(xrange)
	plt = plot(x, Hill.(gate...,x),lw=3, color = :darkgreen, labels = L"Hill(x) = dn + \frac{(up - dn) * K^{n}}{(K^{n} + x^{n})}",
	ylims =[0,8],tickfont = Plots.font("Helvetica", 6), guidefont = Plots.font("Helvetica", 6), legendfontsize = 3, legend = false, fg_legend = :transparent)
	xlabel!("Input");ylabel!("Output")
end
function check_gates_activation_function(xrange,p_unique)
	plt_dynamics = check_δ_range(p_unique, 420)
	plt_gates = [plot_activation(xrange, i) for i in p_unique]
	plt = vcat(plt_dynamics,plt_gates)
	# @show plt
	l = @layout [a b; c d; e f; g h]
	plot(plt...,layout = l)
	# title!("Activation functions")
end

check_gates_activation_function(0:0.01:2, p_unique)

function check_constant_input(p_unique; t_on = 4000, xrange = 0:0.01:2)
	gate_p_set = gate_param_assign(p_unique...)
	# test equilibrium for the example
	u0 =  rand(1:22., length(ode_f1.syms))
	p = 0;tspan = (0.0, 1000.0)
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
	# plot(sol0,vars=[:m_PsrA :m_BetI])
	plt_ecoli_eq = plot(sol0, lw = 2,
	  xlabel = " Time steps", ylabel = "Concentration",
	  color = [:blue :red],
	  vars=[6,7],
	  label =["Q" L"\widebar Q"], # use the  mathemtical gate name
	  title = "Stable States")
	# turn on signal A = 20
	# t_on =3000.0;
	param = [ gate_p_set.g1...,
	       gate_p_set.g2...,
	       gate_p_set.g3...,
	       gate_p_set.g4...,
	       gate_p_set.g5...,
	       gate_p_set.g6...,
	       gate_p_set.g7...,
	       0.025,0.025,20]
	prob1 = ODEProblem(ode_f1, sol0[end], (0.0, t_on), param)
	sol1 = solve(prob1, Tsit5())
	plt_ecoli_const_input = plot(sol1, vars = [6,7], lw = 1.5,
	   xlabel = "time steps (min)", ylabel = "concentration",
	   label =["Q: PsrA" L"$\barQ$: BetI"],
	   title = "Constant input signal",
	   legend = :topright,
	   # ylims = (0.,5.),
	   dpi = 300)
	plt_activations = check_gates_activation_function(xrange,p_unique)
	l = @layout [a b; c]
	return plot(plt_ecoli_eq,plt_ecoli_const_input,plt_activations, layout = l,size=(750,800))
end




pp = check_δ_range(p_unique, 420)

check_constant_input(p_unique, t_on = 3000, xrange = 0:0.01:1)
# check_constant_input(p_unique_cello)



## ----------------------------------------------------------------
# ----check -- the cello produced dataset
# mutable struct Param_E
#     repressor
#     RBS
#     min:: Float64
#     max:: Float64
#     K:: Float64
#     n:: Float64
# end
# para = CSV.File("DEmodels/param_db/para_s4.csv") |> DataFrame
# HlyIIR = Param_E(Matrix(para)[7,:]...)
# LmrA = Param_E(Matrix(para)[10,:]...)
# SrpR = Param_E(Matrix(para)[20,:]...)
# BM3R1 = Param_E(Matrix(para)[6,:]...)
# PhIF = Param_E(Matrix(para)[11,:]...)
# PsrA = Param_E(Matrix(para)[14,:]...)
# BetI = Param_E(Matrix(para)[3,:]...)

# cello_set = [HlyIIR, LmrA, SrpR,BM3R1,PhIF,PsrA,BetI]
# p_unique_cello =[]
# [push!(p_unique_cello,[i.min, i.max, i.K, i.n]) for i in cello_set]
