#Author: Tianchi Chen
## ------ Import package and functions
using DifferentialEquations, ModelingToolkit
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
using Latexify, Random, Base
using CSV, DataFrames, ProgressMeter
# include("functions.jl")# using BlackBoxOptim, LinearAlgebra

## ==== Build multiple counter connectors : 1Bit counter case 📗 =========

hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
 # Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
D = Differential(t)
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

##
# equations(de1)

u0 = rand(7)
tspan = (0.0,1e4)
param = [1.5, 0.002, 0.051, 2.0, 0.025, 0.025, 0]
prob = ODEProblem(ode_f1, u0, tspan, param)
sol0 = solve(prob)
plot(sol0,vars = [m1_HKCI,m1_PhlF])
plot(sol0)



##
h(out, p, t) = (out.=0.0)
const out = zeros(7)
function bit1_delay(du, u, h, p, t)
    m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI, m1_PhlF= u
    up, dn, K, n, γ, ξ, input= p
    h(out, p, t-tau)
    du[1] = ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_PhlF + input)^n))       - γ*out[1]
    du[2] = ξ * (dn + (up - dn) * K^n / (K^n + abs(input + out[1])^n))        - γ*out[2]
    du[3] = ξ * (dn + (up - dn) * K^n / (K^n + abs(out[1] + m1_PhlF)^n))      - γ*m1_CI1
    du[4] = ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_CI1 + out[2])^n))       - γ*m1_PsrA
    du[5] = ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_PsrA)^n))               - γ*m1_BM3RI
    du[6] = ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_BM3RI + m1_PhlF)^n))    - γ*m1_HKCI
    du[7] = ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_HKCI + m1_PsrA)^n))     - γ*m1_PhlF
end
##
# h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(3)
u0 = rand(7)
tspan = (0.0,1e4)
param = [1.5, 0.002, 0.051, 2.0, 0.025, 0.025, 0.0]

prob = ODEProblem(bit1_delay, u0, tspan, param)
sol0 = solve(prob)
plot(sol0,vars = [6,7])
## delayed case without input
tau = 1.1
lags = [tau]
prob = DDEProblem(bit1_delay,u0,h,tspan,param; constant_lags=lags)
alg = MethodOfSteps(Tsit5())
sol = solve(prob,alg)
plot(sol)

##
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

##
u0 = rand(7)
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 250, cycle = 10)
param = [1.1,0.002,0.021,1.5,0.025,0.025,p]
prob0 = ODEProblem(bit1_delay, u0, tspan, param)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
plot(sol, vars = [6, 7])

## delayed case with input
tau = 0.9
lags = [tau]
u0 = [0.69,0.13,0.57,0.43, 0.83, 0.37, 0.53]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 250, cycle = 10)
param = [1.1,0.002,0.021,1.5,0.025,0.025, p]

prob = DDEProblem(bit1_delay,u0,h,tspan,param; constant_lags=lags)
alg = MethodOfSteps(Tsit5())
sol = solve(prob,alg,callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
plot(sol,vars = [6, 7])

##
# using ParameterizedFunctions
#
# de1 = @ode_def delay begin
#     dm1_LexA1 =  ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_PhlF + input)^n))       - γ*m1_LexA1
#     dm1_IcaR  =  ξ * (dn + (up - dn) * K^n / (K^n + abs(input + m1_LexA1)^n))       - γ*m1_IcaR
#     dm1_CI1   =  ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_LexA1 + m1_PhlF)^n))       - γ*m1_CI1
#     dm1_PsrA  =  ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_CI1 + m1_IcaR)^n))       - γ*m1_PsrA
#     dm1_BM3RI =  ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_PsrA)^n))       - γ*m1_BM3RI
#     dm1_HKCI  =  ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_BM3RI + m1_PhlF)^n))    - (γ*m1_HKCI)
#     dm1_PhlF  =  ξ * (dn + (up - dn) * K^n / (K^n + abs(m1_HKCI + m1_PsrA)^n))    - (γ*m1_PhlF)
# end up dn K n γ ξ input
# u0 = rand(7)
# tspan = (0.0,1e4)
# param = [1.5, 0.002, 0.051, 2.0, 0.025, 0.025, 0.0]
# prob = ODEProblem(de1, u0, tspan, param)
# sol0 = solve(prob)
# plot(sol0)
#
#
# u0 = rand(7)
# Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 250, cycle = 10)
# param = [1.1,0.002,0.021,1.5,0.025,0.025,p]
# prob0 = ODEProblem(de1, u0, tspan, param)
# sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
# plot(sol, vars = [6, 7])







## create inset plot
pyplot()
begin
    plt = plot(sol, vars = [6, 7],legend = :topright, label = ["G6" "G7"],
                lw = 1.5, xlabel = "Time", ylabel = "Concentration", dpi = 500)
    scatter!(plt, ts, 1*ones(length(ts)),label="Signal")
    lens!(plt, [19900.0, 20950.0], [0.0, 1.1],  xticks = 19900:200:20950, grid = nothing, title = "Zoomed transition",
    inset = (1, bbox(0.1, 0.2, 0.4, 0.2,:right,:bottom)), subplot = 2)
    savefig(plt, "inset_plot.png")
end
