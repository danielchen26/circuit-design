## Load package
using OrdinaryDiffEq, ModelingToolkit, CSV, Statistics
using DataFramesMeta, DataFrames, LaTeXStrings, Plots
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
        	gate_p_set = gate_param_assign([df[rand_idx[1],:].dn, df[rand_idx[1],:].up, df[rand_idx[1],:].K, df[rand_idx[1],:].n],
                               [df[rand_idx[2],:].dn, df[rand_idx[2],:].up, df[rand_idx[2],:].K, df[rand_idx[2],:].n],
                               [df[rand_idx[3],:].dn, df[rand_idx[3],:].up, df[rand_idx[3],:].K, df[rand_idx[3],:].n],
                               [df[rand_idx[4],:].dn, df[rand_idx[4],:].up, df[rand_idx[4],:].K, df[rand_idx[4],:].n],
                               [df[rand_idx[5],:].dn, df[rand_idx[5],:].up, df[rand_idx[5],:].K, df[rand_idx[5],:].n],
                               [df[rand_idx[6],:].dn, df[rand_idx[6],:].up, df[rand_idx[6],:].K, df[rand_idx[6],:].n],
                               [df[rand_idx[7],:].dn, df[rand_idx[7],:].up, df[rand_idx[7],:].K, df[rand_idx[7],:].n])
							   return gate_p_set, rand_idx
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

@show gate_p_set,rand_idx = gate_p_set_gen(rand(), df; shared = "random")
dump(gate_p_set)
@show rand_idx

## Run one example
"Remember to change the inital control parameter index, here 31 instead of 7.
31 because the 1bit counter contains 31 parameters, the control p hass the last index"
function run_prob_1bit(;init_relax, duration,relax,signal, gate_p_set::gate_param_assign)
    u0 =  rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(index = 31, Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = 10)
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


## set the shared single parameter set index 🔺
idx = 1000
idx_set = [1,2,3,4,5,6,7]
gate_p_set = gate_p_set_gen(idx_set, dff, shared="gaussian")
sol, ts = run_prob_1bit(;init_relax = 5000., duration=dff[Int64(median(idx_set)),:].δ, relax=5000., signal=20., gate_p_set);
plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"])



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
##


## test random selected 3 gate from para_s4 dataset to see if counter works.
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
# convert.(Int64, rand_idx) = [18, 9, 17, 18, 18, 18, 18]    340<δ<430
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
specific_idx = [18, 9, 17, 18, 18, 18, 18]
p_set_specific = []
[push!(p_set_specific,[df[specific_idx[i],:].dn, df[specific_idx[i],:].up, df[specific_idx[i],:].K, df[specific_idx[i],:].n]) for i in 1:7]
gate_p_set_specific = gate_param_assign(p_set_specific...)
sol, ts = run_prob_1bit(;init_relax = 5000., duration= 340, relax=5000., signal=20., gate_p_set = gate_p_set_specific);
plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = L"\delta=\ %$dt", ylims = (0.,3))
display(plt)
# check cost function
costtot = cost_bit1(sol, ts, up)

## test random selected 7 gate from para_s4 dataset to see if counter works.
rand_idx_δ_set = []
gates_p_set = []

for i in 1:30
	gate_p_set,rand_idx = gate_p_set_gen(rand(), df; shared = "random")
	dump(gate_p_set)
	all_7_up_set = [getfield(gate_p_set,name)[2] for name in fieldnames(gate_param_assign)]
	up = mean(all_7_up_set)
	for dt in 200:5:500
		sol, ts = run_prob_1bit(;init_relax = 5000., duration= dt, relax=5000., signal=20., gate_p_set);
		costtot = cost_bit1(sol, ts, up)
		@show costtot, dt
		costtot == 0 ? push!(rand_idx_δ_set, [rand_idx,dt]) && push!(gates_p_set,gate_p_set) : nothing
		plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],label =["Q" L"\overline{Q}"], title = "$dt,$rand_idx")
		display(plt)
	end
end


dt_set










## test regular 1bit counter
# function run_prob_1bit(;init_relax, duration,relax,signal,K,n,up)
#     u0 =  rand(1:22., length(ode_f1.syms))
#     Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = 20)
#     param = [up,0.002,K,n,0.025,0.025,p]
#     prob0 = ODEProblem(ode_f1, u0, tspan, param)
#     sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
#     return sol, ts
# end
# df_1bit = DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[], up = Float64[])
# tt = []
# @time @showprogress for K = 0.01: 0.05:1., n = 1.:0.2:10., δ = 10:5:500, A = 20., up = 1:1.:10 #@showprogress
#     # for K = 0.011: 0.05:1., n = 1.5:0.5:10., δ = 250:10:350, A = 20., up = 1:0.5:10 # try broader range
#     # K = 0.011: 0.02:0.1, n = 1.5:0.1:3., δ = 250:5:350, A = 20., up = 1:0.1:3, # original
#     # solve DiffEq given parameters
#     sol, ts = run_prob_1bit(;init_relax = 2500., duration=δ, relax=2000., signal=A, K=K, n=n, up = up)
#     # Get cost for this parameters set
#     costtot = cost_bit1(sol, ts, up)
#     println("K:$K n:$n, δ:$δ, A:$A, up:$up\n")
#     @show costtot
#     costtot == 0 ? push!(df_1bit, [K, n, δ, A, up]) : nothing
#
#     # count example
#     push!(tt,1.)
#     @show sum(tt)
#     # if sum(tt) >= 800.
#     #     break
#     # end
# end




## Sampling two different random gates from library, 5 identical
