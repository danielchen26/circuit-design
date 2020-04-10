## Import Functions
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
## ==== Build multiple counter connectors : 3 Bits counter case 📗 =========
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
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t) g1(t) g2(t) g3(t) m2_LexA1(t) m2_IcaR(t) m2_CI1(t) m2_PsrA(t) m2_BM3RI(t) m2_HKCI(t) m2_PhlF(t) g21(t) g22(t) g23(t) g24(t) g25(t) g26(t) m3_LexA1(t) m3_IcaR(t) m3_CI1(t) m3_PsrA(t) m3_BM3RI(t) m3_HKCI(t) m3_PhlF(t)
@derivatives D'~t
##
eqs = [
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
##
de = ODESystem(eqs, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF, g1, g2, g3, m2_LexA1, m2_IcaR, m2_CI1, m2_PsrA,m2_BM3RI, m2_HKCI, m2_PhlF, g21, g22, g23, g24, g25, g26, m3_LexA1, m3_IcaR, m3_CI1, m3_PsrA, m3_BM3RI, m3_HKCI, m3_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f = ODEFunction(de)
## Solve DifferentialEquations
u0 =  Float64[i for i in rand(1:22, 30)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 270, cycle = 20)
param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
prob0 = ODEProblem(ode_f, u0, tspan, param)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
## Plots
# up = 1.5
# Ylim = 3*up
# # plt_opt = [legend = :topright, legendfontsize = 7, ylims =(0.,Ylim)]
# C_plt = plot(sol, vars = [:g3], legend = false, ylims =(0.,Ylim))
# C2_plt = plot(sol, vars = [:g26], legend = false, ylims =(0.,Ylim))
# B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = false, ylims =(0.,Ylim))
# B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = false, ylims =(0.,Ylim))
# B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = false, ylims =(0.,Ylim))
# Bit2_plt = plot(C_plt,C2_plt,B1_plt,B2_plt,B3_plt,layout = (5,1),size=(600,500))
## optimization for parameters searching

## Wrap the above problem to a fucntion
function run_prob_3bits(;init_relax,duration,relax,signal,K,n,up,cycle)
    u0 =  Float64[i for i in rand(1:22, 30)]
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = cycle)
    param = [up,0.002,K,n,0.025,0.025,p]
    prob0 = ODEProblem(ode_f, u0, tspan, param)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
    return sol, ts
end

sol, ts = run_prob_3bits(;init_relax = 20000., duration=270.,relax=20000.,signal=20.,K=0.081,n=2.81, up= 1.5, cycle=20)
# sol, ts = run_prob_3bits(;duration=285.,relax=20000.,signal=10.,K=0.091,n=2.9)
##
# Ylim = 3*up
# δ = 270.
# sol,ts= run_prob_3bits(duration = δ, relax = 20000., signal = 20.)
# plot(sol, vars = [:m_HKCI, :m_PhlF, :m2_HKCI, :m2_PhlF, :m3_HKCI, :m3_PhlF ], legend = false, layout=(3,2), ylims =(0.,Ylim), title=["Bit1" "Bit1" "Bit2" "Bit2" "Bit3" "Bit3"], title_location=:left)


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
ts_ID, C1, C2  = Carrying_bit(sol, ts)

function MB_T0(ts_ID, C1, C2, up)
    thred =up/2
    # Identify the T0
    # C1_B = C1 .> thred; C2_B = C2 .> thred
    T0 = []
    ts_cycle = collect(1:2:length(ts_ID))
    for i = 1:length(C1)
        if C1[i] > thred && C2[i] > thred
            # @show i , C1[i], C2[i], thred
            push!(T0,  ts_cycle[i], ts_cycle[i] + 1)
            break
        else
            continue
        end
    end

    return T0
end
# Identify the T0 for multibit counter
T0 = MB_T0(ts_ID, C1, C2, up) # if something goes wrong, this may give me 0

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


function osci(G6, up)
    thred = up/2
    G6b = G6 .> thred
    G6b1 = [G6b[i] for i = 1:2:length(G6b)]
    G6b2 = [G6b[i] for i = 2:2:length(G6b)]
    sum(G6b1 .* G6b2)
end


function cost_bit3(sol, ts, up)
    thred = up/2
    # set T0 for multibit counter
    ts_ID, C1, C2  = Carrying_bit(sol, ts)
    T0 = MB_T0(ts_ID, C1, C2, up)
    @show T0
    if length(T0) > 0
        G6,  G7 = Switch_cost(sol, [6,7], Δ0, ts, T0)
        G16, G17 = Switch_cost(sol, [16,17], Δ0, ts, T0)
        G29, G30 = Switch_cost(sol, [29,30], Δ0, ts, T0)
        # @show G6
        @show osci(G6, up), osci(G7, up)
        if osci(G6, up) > 0 || osci(G7, up) > 0
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
costtot = cost_bit3(sol, ts, up)
##

## Parameters Searching for  3bits counter
# ======== Sampling Parameters ===================
# Varying K, n, δ
df3 = DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[], up = Float64[])
tt = []
@time @showprogress for K = 0.001: 0.02:0.1, n = 1.5:0.1:4., δ = 250:5:350, A = 20., up = 1:0.1:3 #@showprogress
    # solve DiffEq given parameters
    sol, ts = run_prob_3bits(;init_relax = 20000.,duration=δ, relax=20000., signal=A, K=K, n=n, up = up, cycle = 20)
    # Get cost for this parameters set
    costtot = cost_bit3(sol, ts, up)
    println("K:$K n:$n, δ:$δ, A:$A, up:$up\n")
    @show costtot
    costtot == 8 ? push!(df3, [K, n, δ, A, up]) : nothing

    # count example
    push!(tt,1.)
    @show sum(tt)
    if sum(tt) >= 8.
        break
    end
end
##
df3
CSV.write("3Bits_DB.csv", df3)
