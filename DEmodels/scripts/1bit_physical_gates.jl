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
# the best circuit for 1 bit that has the largest cv among 7 circuits.
p_unique = Any[[0.03, 2.8, 0.21, 2.4],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]



function check_δ_range(p_unique, gate_out_names, δ_set; Δ0 = 1500, Δ = 1500, cycle = 8, legendfontsize = 3, signal = 10, Ylim = 15)
	gate_p_unique = gate_param_assign(p_unique...)
	if length(δ_set) ==1
		sol, ts = run_prob_1bit(;init_relax = Δ0, duration = δ_set, relax = Δ, signal=signal, gate_p_set = gate_p_unique, cycle = cycle);
		@show costtot = cost_bit1(sol, ts)
		plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],lw = 1.5,
          		   xlabel = "time steps (min)", ylabel = "concentration",
				#    label =["Q" L"\overline{Q}"],
				#    label =["Q: PsrA" L"$\barQ$: BetI"], # note that these two specific gate need to match the p_unique order, double check !!
				   label = ["Q: "*gate_out_names[1] L"$\barQ$: "*gate_out_names[2]],
				   legendfontsize = legendfontsize, dpi = 500,
				   title = "Circuit dynamics (e.g., signal duration: $δ_set min)")
				#    title = L"Circuit\ dynamics\ (e.g., signal\ duration\ \delta :\ %$δ_set)")
		ts1 = vcat(ts, sol.t[end])
		plot!(plt,ts1, push!(repeat([0,signal],cycle+1),0),
      			linetype=:steppre, linestyle = :dot,
      			color = :purple,
				label = "Input signal",
				ylims =(0.,Ylim) )
		display(plt)
		return plt, costtot
	else
		nothing
	end

	cost_set =[]
	for δ in δ_set
		sol, ts = run_prob_1bit(;init_relax = Δ0, duration = δ, relax = Δ, signal=20., gate_p_set = gate_p_unique, cycle = cycle);
		# g67min =  maximum([minimum(sol[6,Int64(round(length(sol)/2)):end]),minimum(sol[7,Int64(round(length(sol)/2)):end])]);
		# g67max = minimum([maximum(sol[6,Int64(round(length(sol)/2)):end]),maximum(sol[7,Int64(round(length(sol)/2)):end])]);
		# up = (g67min + g67max)/2
		# @show up
		println("δ values is:",δ)
		costtot = cost_bit1(sol, ts)
		push!(cost_set,costtot)
		plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],
				   label =["Q" L"\overline{Q}"],
				   legendfontsize = legendfontsize, 
				   title = L"Circuit\ dynamics (e.g., signal\ duration\ \delta :\ %$δ")
		display(plt)
		
		# if length(δ_set) ==1
		# 	return plt, costtot
		# else
		# 	nothing
		# end
	end
	return cost_set
end

Hill(dn, up, K, n, x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
function plot_activation(xrange,gate)
	x = collect(xrange)
	plt = plot(x, Hill.(gate...,x),lw=3, 
			   color = :darkgreen, 
			   scale=:log10,
			   labels = L"Hill(x) = dn + \frac{(up - dn) * K^{n}}{(K^{n} + x^{n})}",
	ylims =[0,8],tickfont = Plots.font("Helvetica", 6), guidefont = Plots.font("Helvetica", 6), legendfontsize = 3, legend = false, fg_legend = :transparent)
	xlabel!("Input");ylabel!("Output")
end
function check_gates_activation_function(xrange,p_unique)
	plt_dynamics, costtot = check_δ_range(p_unique, 420)
	plt_gates = [plot_activation(xrange, i) for i in p_unique]
	plt = vcat(plt_dynamics,plt_gates)
	# @show plt
	l = @layout [a b; c d; e f; g h]
	return plot(plt...,layout = l), costtot
	# title!("Activation functions")
end
# check_gates_activation_function(0:0.01:2, p_unique)

function check_all(p_unique; t_on = 4000, xrange = 0:0.01:2)
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
	plt_activations,costtot = check_gates_activation_function(xrange,p_unique)
	l = @layout [a b; c]
	plt_check_all = plot(plt_ecoli_eq,plt_ecoli_const_input,plt_activations, layout = l,size=(800,800))
	return plt_check_all, costtot
end

plot_activation(0:0.01:2,p_unique[1])

# ---------- test some plots for activation functions ↓
x = 1e-3:0.01:100
using GRUtils
gr() # the following plot need gr() backend.
# inspectdr()
plot(x, Hill.(p_unique[2]...,x),lw=10, 
			   color = colorant"#a0c54d", 
			   xaxis=:log, yaxis=:log,
			   thickness_scaling = 1.3,
			#    bordercolor="white",
			#    tickfontsize = [],
			# ticklabels = [],

			   framestyle = :box,
			   grid=false,
			#    foreground_color_tick = :green, 
			#    foreground_color_minortick=:red, 
			   minorgrid=true,
			   minorticks = 10,
			   labels = L"Hill(x) = dn + \frac{(up - dn) * K^{n}}{(K^{n} + x^{n})}",
	ylims =[1e-2,8],tickfont = Plots.font("Helvetica", 6), guidefont = Plots.font("Helvetica", 6), legendfontsize = 3, legend = false, fg_legend = :transparent)
	xlabel!("Input");ylabel!("Output")

# using Colors
# colorant"#a0c54d"
# ---------- test some plots for activation functions ⤉

# below is to save one of the best circuit among 7 physical_gate_cases and used in figure 5 in paper.
plt_1_7_best, cost_set = check_δ_range(p_unique, 450, legendfontsize = 8, Ylim = 13)
# savefig(plt_1_7_best,"./DEmodels/scripts/Paper_plots/physical_gate_cases/best1bit_from7(4th).png")

plt_check_all, costtot = check_all(p_unique, t_on = 3000, xrange = 0:0.01:1)
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

# cello_set = [HlyIIR,LmrA,SrpR,BM3R1,PhIF,PsrA,BetI]
# p_unique_cello =[]
# [push!(p_unique_cello,[i.min, i.max, i.K, i.n]) for i in cello_set]


# ===== Table:replace ======
# if to replace 1 cello gate
# │ 1   │ AmeR      │ F1     │ 0.2     │ 3.8     │ 0.09    │ 1.4     │
# │ 2   │ AmtR      │ A1     │ 0.06    │ 3.8     │ 0.07    │ 1.6     │
# │ 8   │ IcaRA     │ I1     │ 0.08    │ 2.2     │ 0.1     │ 1.4     │
# │ 9   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 15  │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 16  │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │

replace_gates = [[0.2  ,3.8  ,0.09 ,1.4 ],[0.06 ,3.8  ,0.07 ,1.6 ],[0.08 ,2.2  ,0.1  ,1.4 ],[0.07 ,4.3  ,0.05 ,1.7 ],[0.01 ,2.4  ,0.05 ,2.7 ],[0.03 ,2.8  ,0.21 ,2.4 ]]

# change 1 cello gate
p_unique_cello = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
pp = check_δ_range(p_unique, 200:5:500)
check_all(p_unique_cello1, t_on = 3000, xrange = 0:0.01:1)
@show DataFrame(p_unique_cello)

# test for each row in Table:replace above, if the replacement of one row in cello circuit is still able to produce a feasible counter.
modi_1gate = []
for i ∈ replace_gates
	for j in 1:7
		# @show i
		modi = copy(p_unique_cello)
		modi[j] = i
		@show DataFrame(modi)
		plt_check_all, costtot = check_all(modi, t_on = 3000, xrange = 0:0.01:1)
		println("Cost function:", costtot)
		display(plt_check_all)
		if costtot == 0
			push!(modi_1gate, DataFrame(modi))
			savefig(plt_check_all, "./DEmodels/scripts/Paper_plots/physical_gate_cases/cello_replace_"*"$i"*"_"*"$j"*".png")
		end
		# @show DataFrame(modi)
		# println("\n")
	end
end




## Genereate 7 circuits' δ ranges and choose the best δ with largest range.
###########################################
function load_experiment_gates_circuits()
	# 1. cello  gates in order: => [HlyIIR, LmrA, SrpR, BM3R1, PhIF, PsrA,BetI]
	p_unique1 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
	circuit1_name = ["H1HlyIIR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
	# 2. cello_replace_[0.2, 3.8, 0.09, 1.4]_2   🍏 LmrA ==> AmeR
	p_unique2 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 3.8, 0.09, 1.4],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
	circuit2_name = ["H1HlyIIR", "F1AmeR", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
	# 3. cello_replace_[0.01, 2.4, 0.05, 2.7]_5  🍏 PhIF ==> QacR
	p_unique3 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 2.4, 0.05, 2.7],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
	circuit3_name = ["H1HlyIIR", "N1LmrA", "S4SrpR", "B3BM3R1", "Q1QacR", "R1PsrA","E1BetI"]
	# 4. 🔴 cello_replace_[0.03, 2.8, 0.21, 2.4]_1  🍏 HlyIIR ==> QacR
	p_unique4 = Any[[0.03, 2.8, 0.21, 2.4],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
	circuit4_name = ["Q2QacR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
	# 5. cello_replace_[0.03, 2.8, 0.21, 2.4]_7  🍏 BetI ==> QacR
	p_unique5 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.03, 2.8, 0.21, 2.4]]
	circuit5_name = ["H1HlyIIR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","Q2QacR"]
	# 6. cello_replace_[0.06, 3.8, 0.07, 1.6]_1  🍏 HlyIIR ==> AmtR
	p_unique6 = Any[[0.06, 3.8, 0.07, 1.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
	circuit6_name = ["A1AmtR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
	# 7. cello_replace_[0.07, 4.3, 0.05, 1.7]_1  🍏 HlyIIR ==> LitR
	p_unique7 = Any[[0.07, 4.3, 0.05, 1.7],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
	circuit7_name = ["L1LitR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]

	p_unique_set = [p_unique1,p_unique2,p_unique3,p_unique4,p_unique5, p_unique6,p_unique7]
	circuit_name = [circuit1_name, circuit2_name, circuit3_name, circuit4_name, circuit5_name, circuit6_name, circuit7_name]
	return p_unique_set, circuit_name
end 

# load each circuit names and hill parameters for the 7 feasible counters with unique gates
p_unique_set, circuit_name = load_experiment_gates_circuits()

# Generate all cost functions values for each circuit shown above (7 in total)
all_cost = []
for p_unique in p_unique_set
	@show DataFrame(p_unique)
	cost_set_i = check_δ_range(p_unique, 200:5:600)
	push!(all_cost,cost_set_i)
end 
df_7circuits = DataFrame(all_cost)
df_7circuits = df_7circuits .==0
df_7circuits.δrange = 200:5:600
showall(df_7circuits)
CSV.write("./df_7circuits_cost.csv", df_7circuits)

new_df = CSV.File("./df_7circuits_cost.csv") |> DataFrame
new_names = ["Circuit 1","Circuit 2","Circuit 3","Circuit 4","Circuit 5","Circuit 6","Circuit 7","δrange"]
rename!(new_df, new_names)

test = stack(new_df,1:7)
new_test = test[test.value .==1,:]
new_test |>
	@vlplot(
		width=500,
		height=200,
		mark={:boxplot},
		x={"variable:o", axis={title="Circuit NO."} },
		y={:δrange, scale = {domain = (290, 570)}},
		color={"value", scale={domain=["variable"],range=["#393b79"]}, legend = nothing}
		) |> save("7circuit_δ_boxplot.svg")







# Generate dynamics plots for 7 circuits(with unique experimental gates) with δ = 400
for i in 1:7
	p_unique = p_unique_set[i]
	gate_out_names = circuit_name[i][[6,7]]
	println("The circuit gates parameters: ", p_unique)
	println("The gates for this circuit are: ", circuit_name[i])
	println("The last output gates are: ", gate_out_names)
	δ = 400
	plt, costtot = check_δ_range(p_unique, gate_out_names, δ,legendfontsize = 10)
	savefig(plt,pwd()*"/DEmodels/scripts/Paper_plots/physical_gate_cases/circuit_$i/dynamics_δ = $δ"*".png")
	println("Cost: ", costtot)
	println("\n")
end


# ======= Generate the following table, need `df_7circuits` from above ====

# │ Row │ Circuit │ δ_min │ δ_max │ δ_range │ δ_cv     │
# │     │ Int64   │ Int64 │ Int64 │ Int64   │ Float64  │
# ├─────┼─────────┼───────┼───────┼─────────┼──────────┤
# │ 1   │ 1       │ 375   │ 460   │ 85      │ 0.184783 │
# │ 2   │ 2       │ 360   │ 470   │ 110     │ 0.234043 │
# │ 3   │ 3       │ 325   │ 505   │ 180     │ 0.356436 │
# │ 4   │ 4       │ 300   │ 550   │ 250     │ 0.454545 │
# │ 5   │ 5       │ 385   │ 455   │ 70      │ 0.153846 │
# │ 6   │ 6       │ 320   │ 515   │ 195     │ 0.378641 │
# │ 7   │ 7       │ 345   │ 490   │ 145     │ 0.295918 │

range_set =[]
for i in 1:7
	range = df_7circuits[df_7circuits[:,i].==0,:].δrange[[1,end]]
	push!(range_set,range)
end
range_set
circuits_range = DataFrame(Circuit = 1:7, δ_min = [ i[1] for i in range_set], δ_max = [ i[2] for i in range_set])
circuits_range.δ_range = circuits_range.δ_max - circuits_range.δ_min
circuits_range.δ_cv = circuits_range.δ_range ./circuits_range.δ_max
circuits_range


using VegaLite
circuits_range |>
        @vlplot(
            width=500,
            height=200,
            mark={:boxplot, extent="min-max"},
            x="Circuit:o",
            y={:δ_range, axis={title="Duration δ"}})















## Now test the case | work here to select more feasible circuit other than cello one.
################################################################
gate_p_set,df_unique = gate_p_set_gen(rand(), df; shared = "random", rand_num = 7)

p_unique = [];[push!(p_unique,[i.dn, i.up, i.K, i.n]) for i in eachrow(df_unique[:, [:dn,:up,:K, :n]])]
@show p_unique


@show df_unique
@show filter(row -> row.repressor ∉ df_unique.repressor, df)

replace_gates_df = filter(row -> row.repressor ∉ df_unique.repressor, df)
replace_gates_df[:,[:dn,:up,:K, :n]]

modi_1gate = []
for i ∈ replace_gates_df
	for j in 1:7
		# @show i
		modi = copy(p_unique)
		modi[j] = i
		@show DataFrame(modi)
		plt_check_all, costtot = check_all(modi, t_on = 3000, xrange = 0:0.01:1)
		println("Cost function:", costtot)
		display(plt_check_all)
		if costtot == 0
			push!(modi_1gate, DataFrame(modi))
			savefig(plt_check_all, "./DEmodels/scripts/Paper_plots/physical_gate_cases/test1_replace_"*"$i"*"_"*"$j"*".png")
		end
	end
end











## jacobian
