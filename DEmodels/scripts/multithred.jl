
## ------ Import package and functions
using ModelingToolkit, DifferentialEquations, Plots
using Optim
using Printf, DataFrames, CSV
# using ProgressMeter,
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

## Test parallel computation
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
### Define a differential equation system
@parameters t up dn K n γ ξ p
@variables g1(t) g2(t) g3(t) g4(t) g5(t) g6(t) g7(t)
@derivatives D'~t
##
eqs1 = [
    # Bit 1 =================
    D(g1) ~ ξ * hill(g7 + p)        - deg(g1),
    D(g2 ) ~ ξ * hill(g1 + p)       - deg(g2),
    D(g3  ) ~ ξ * hill(g1 + g7) - deg(g3),
    D(g4 ) ~ ξ * hill(g2 + g3)   - deg(g4),
    D(g5) ~ ξ * hill(g4)            - deg(g5),
    D(g6 ) ~ ξ * hill(g5 + g7) - deg(g6),
    D(g7 ) ~ ξ * hill(g4 + g6)  - deg(g7)]
##
de1 = ODESystem(eqs1, t, [g1, g2, g3, g4, g5, g6,g7], [up,dn,K,n,γ,ξ,p])
ode_f1 = ODEFunction(de1)

## Run the 1bit counter problem
# # ==== Randomized Initials equlibrations =====
u0 = rand(1:22., length(ode_f1.syms))
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, pp = init_control(Δ0 = 1000., Δ = 1000., δ = 270., A = 20., cycle = 10)
param = [1.5,0.002,0.081,2.81,0.025,0.025,pp]
prob0 = ODEProblem(ode_f1, u0, tspan, param)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts)
plot(sol, vars=[:g6,:g7])
# ## parallel working
#
# ptest = [[up, dn, k, n] for up = 1.5 for dn = 0.002 for k = 0.011:0.01:0.081 for n = 1:0.1:2.9]
# ptest
# function prob_func(prob,i,repeat)
#     # param = [1.5,0.002,0.081,nn[i],0.025,0.025,pp]
#     param = [ptest[i]..., 0.025,0.025,pp]
#     # println("parameters", param)
#     remake(prob; p = param)
# end
# ensemble_prob = EnsembleProblem(prob0,prob_func=prob_func)
# @time sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=160, callback = cb, tstops = ts)
#
# ##
#
# for i = 1:160
#     param =  sim[i].prob.p
#     display(plot(sim[i], vars = [:g6, :g7],title = "Param: $param") )
# end
## complete parellel prob define
##
function complete_parellel_prob(;δ , A )
    # define single problem
    u0 = rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, pp = init_control(Δ0 = 5000., Δ = 5000., δ = δ, A = A, cycle = 20)
    param = [1.5,0.002,0.081,2.81,0.025,0.025,pp]
    prob0 = ODEProblem(ode_f1, u0, tspan, param)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts)

    # multithreding problem
    ptest = [[up, dn, k, n] for up = 1.5:0.2:3. for dn = 0.002 for k = 0.011:0.01:0.081 for n = 1:0.1:3]
    ptest
    function prob_func(prob,i,repeat)
        # param = [1.5,0.002,0.081,nn[i],0.025,0.025,pp]
        param = [ptest[i]..., 0.025,0.025,pp]
        # println("parameters", param)
        remake(prob; p = param)
    end
    ensemble_prob = EnsembleProblem(prob0,prob_func=prob_func)
    @time sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=1344, callback = cb, tstops = ts)
end


##
# using BenchmarkTools
# @btime complete_parellel_prob(;δ=270. , A = 20.)


##
sim_set = []
for δ = 270:5:280, A = 10:10:30
    @time sim = complete_parellel_prob(;δ =δ, A = A )
    push!(sim_set,sim)
    GC.gc()
end

sim_set

sim_set[1][1].prob.p



##
# multithreding problem
ptest = [[up, dn, k, n, init_control(Δ0 = 5000., Δ = 5000., δ = δ, A = A, cycle = 20)[9], init_control(Δ0 = 5000., Δ = 5000., δ = δ, A = A, cycle = 20)[10] ] for up = 1.5:0.2:3. for dn = 0.002 for k = 0.011:0.01:0.081 for n = 1:0.1:3 for δ = 270:5:280 for A = 10:10:30 ]
ptest

function prob_func(prob,i,repeat)
    # param = [1.5,0.002,0.081,nn[i],0.025,0.025,pp]
    param = [ptest[i][1:4]..., 0.025,0.025,pp]
    # println("parameters", param)
    cb = cb_set[i]
    ts = ts_set[i]
    remake(prob; p = param, callback = ptest[i][6], tstops = ptest[i][5])
end
ensemble_prob = EnsembleProblem(prob0,prob_func=prob_func)
@time sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=12096)
##







@code_warntype complete_parellel_prob(;δ =δ, A = A )











up=1.5
solss =[]
@btime for  n = range(1, stop=4, length=10)
    u0 = rand(1:22., length(ode_f1.syms))
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, pp = init_control(Δ0 = 1000., Δ = 1000., δ = 270., A = 20., cycle = 10)
    param = [1.5,0.002,0.081,n,0.025,0.025,pp]
    prob0 = ODEProblem(ode_f1, u0, tspan, param)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts)
    push!(solss,sol)
    #
    # costtot = cost_bit1(sol, ts, up)
    # println("K:$K n:$n, δ:$δ, A:$A, up:$up\n")
end

##

using BenchmarkTools

plot(sim)

v = sim.u[10].u
rslt = hcat(v...)'
plot(sim.u[10].t,rslt[:,[6,7]])
