## Import Functions
using ModelingToolkit, DifferentialEquations, Plots;pyplot()
using Optim
using Printf, DataFrames, CSV, LaTeXStrings
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


## ======= Build multiple counter connectors : 2 Bits counter case 📗 =========
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
### Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t) g1(t) g2(t) g3(t) m2_LexA1(t) m2_IcaR(t) m2_CI1(t) m2_PsrA(t) m2_BM3RI(t) m2_HKCI(t) m2_PhlF(t)
@derivatives D'~t
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

de2 = ODESystem(eqs2, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF, g1, g2, g3, m2_LexA1, m2_IcaR, m2_CI1, m2_PsrA,m2_BM3RI, m2_HKCI, m2_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f2 = ODEFunction(de2)
# ##
# u0 = Float64[i for i in rand(1:22, 17)]
# Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 270, cycle = 10)
# prob0 = ODEProblem(ode_f2, u0, tspan, p)
# sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)

# ## Run the 2bit counter problem
# # # ==== Randomized Initials equlibrations =====
# u0 = rand(1:22., length(ode_f2.syms))
# p = 0.0
# tspan = (0.0, 3000.0)
# param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
# prob0 = ODEProblem(ode_f2, u0, tspan, param)
# sol0 = solve(prob0, Tsit5())
# plot(sol0, lw = 2, ylims = (0, 5))
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


## Paper plot : Figure 10.
# ----- pyplot setting ------
# Plots.scalefontsizes(0.9) # changing overall fontsize
# 2 bits plot
function Bit2_plot(sol,up)
    Ylim = 3*up
    C_plt = plot(sol, vars = [:g3], legend = :topright, 
                label = "Connector type 1",
                legendfontsize	= 7,
                legendtitle = "Connector module",
                ylims =(0.,Ylim))
    B1_plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],
                label = ["Q1: output" L"$\bar{Q}$1: variable to feed back"],
                legend = :topright, legendfontsize	= 7, 
                legendtitle = "The first bit counter",
                ylims =(0.,Ylim))
    scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal",ylims =(0.,Ylim))
    B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],
                label = ["Q2: output" L"$\bar{Q}$2: variable to feed back"],
                legend = :topright, legendfontsize	= 7, 
                legendtitle = "The second bit counter",
                ylims =(0.,Ylim))
    Bit2_plt = plot(
                    C_plt,B1_plt,B2_plt,
                    # B1_plt,C_plt,B2_plt,
                    lw = 2.,layout = (3,1),xtickfontsize=15,ytickfontsize=15,
                    size=(1000,800),legend = :topright, legendfont = font("Times new roman", 12),
                    xlabel = "Time", ylabel = "Concentration",dpi = 400,
                    # title = ["Carrying Bit" "Bit 1" "Bit 2"]
                    # title = ["Bit 1" "Carrying Bit" "Bit 2"]
                    )
    # display(Bit2_plt)
end

K = 0.011; n = 1.5; δ = 300; A = 20.; up = 3.; Δ0 = 5000.
sol, ts = run_prob_2bits(;init_relax = Δ0, duration=δ, relax=Δ0, signal=A, K=K, n=n, up = up, cycle = 7)
plt = Bit2_plot(sol,2.2)
savefig(plt,"./DEmodels/scripts/Paper_plots/2bitdynamics_short")
#
# anim = @animate for up = 2:0.1:3
#     # up = 1.5
#     sol, ts = run_prob_2bits(;init_relax = 2000., duration=270.,relax=2000.,signal=20,K=0.031,n=1.5, up= up,cycle = 20)
#     plt = Bit2_plot(sol,up)
# end
# gif(anim, "bit2_fps15.gif", fps = 1)
#


## Building cost function for param
# up = 1.8
# sol, ts = run_prob_2bits(;init_relax = 5000., duration=275.,relax=5000.,signal=20,K=0.011,n=1.5, up= 1.8,cycle = 20)

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
# ts_ID, C1 = Carrying_bit(sol, ts)

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

function osci(G6, up)
    thred = up/2
    G6b = G6 .> thred
    G6b1 = [G6b[i] for i = 1:2:length(G6b)]
    G6b2 = [G6b[i] for i = 2:2:length(G6b)]
    # @show G6b1, G6b2, G6b1 .+ G6b2
    return sum(G6b1 .+ G6b2), sum(G6b1 .* G6b2) # comp should be 4, and sum should be 0
end


function cost_bit2(sol, ts, Δ0, up)
    thred = up/2
    # set T0 for multibit counter
    ts_ID, C1 = Carrying_bit(sol, ts)
    T0 = MB_T0(ts_ID, C1)
    @show T0
    if length(T0) > 0
        G6,  G7 =  Switch_cost(sol, Δ0, ts, T0; BNO = "B1")
        G16, G17 = Switch_cost(sol, Δ0, ts, T0; BNO = "B2")
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

function plot_2bit(sol, up, param)
    plt0 = plot(sol, vars = [:g3])
    plt1 = plot(sol, vars =[6,7])
    plt2 = plot(sol, vars =[16,17])
    plt = plot(plt0,plt1,plt2, layout = (3,1), ylims =(0.,3*up), title = ["Param: $param" "B1" "B2"])
    display(plt)
end

## Parameters Searching for  2bits counter
# ======== Sampling Parameters ===================
using ProgressMeter
# Varying K, n, δ
df2 = DataFrame(K = Float64[], n = Float64[], δ =  Float64[], A =  Float64[], up = Float64[])
tt = []
@time @showprogress for K = 0.011: 0.02:0.1, n = 1.5:0.1:3., δ = 250:5:350, A = 20., up = 1:0.1:3, Δ0 = 5000. #@showprogress
# @time @showprogress for K = 0.011, n = 1.5, δ = 270, A = 20., up = 1.3:0.1:2., Δ0 = 5000.  # this line is for test
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

# #unused#::Core.Compiler.Const(var"#complete_parellel_prob##kw"(), false)
# @_2::NamedTuple{(:δ, :A),Tuple{Float64,Float64}}
# @_3::Core.Compiler.Const(complete_parellel_prob, false)
# δ::Float64
# A::Float64
# @_6::Float64
# @_7::Float64
#
# Body::EnsembleSolution{_A,_B,_C} where _C where _B where _A
# 1 ─ %1  = Base.haskey(@_2, :δ)::Core.Compiler.Const(true, false)
# │         %1
# │         (@_6 = Base.getindex(@_2, :δ))
# └──       goto #3
# 2 ─       Core.Compiler.Const(:(Core.UndefKeywordError(:δ)), false)
# └──       Core.Compiler.Const(:(@_6 = Core.throw(%5)), false)
# 3 ┄       (δ = @_6)
# │   %8  = Base.haskey(@_2, :A)::Core.Compiler.Const(true, false)
# │         %8
# │         (@_7 = Base.getindex(@_2, :A))
# └──       goto #5
# 4 ─       Core.Compiler.Const(:(Core.UndefKeywordError(:A)), false)
# └──       Core.Compiler.Const(:(@_7 = Core.throw(%12)), false)
# 5 ┄       (A = @_7)
# │   %15 = (:δ, :A)::Core.Compiler.Const((:δ, :A), false)
# │   %16 = Core.apply_type(Core.NamedTuple, %15)::Core.Compiler.Const(NamedTuple{(:δ, :A),T} where T<:Tuple, false)
# │   %17 = Base.structdiff(@_2, %16)::Core.Compiler.Const(NamedTuple(), false)
# │   %18 = Base.pairs(%17)::Core.Compiler.Const(Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}}(), false)
# │   %19 = Base.isempty(%18)::Core.Compiler.Const(true, false)
# │         %19
# └──       goto #7
# 6 ─       Core.Compiler.Const(:(Base.kwerr(@_2, @_3)), false)
# 7 ┄ %23 = Main.:(var"#complete_parellel_prob#22")(δ, A, @_3)::EnsembleSolution{_A,_B,_C} where _C where _B where _A
# └──       return %23
# 13.782367 seconds (30.98 M allocations: 3.884 GiB, 26.72% gc time)
# 9.230517 seconds (30.19 M allocations: 3.874 GiB, 31.18% gc time)
# 11.361198 seconds (31.41 M allocations: 4.023 GiB, 40.28% gc time)
# 9.779346 seconds (30.75 M allocations: 3.942 GiB, 44.75% gc time)
# 9.785 s (30771750 allocations: 3.95 GiB)
# 19.741373 seconds (30.14 M allocations: 3.868 GiB, 61.66% gc time)
# 13.363160 seconds (31.17 M allocations: 3.993 GiB, 40.77% gc time)
# 10.703867 seconds (29.75 M allocations: 3.820 GiB, 41.73% gc time)
# 8.197702 seconds (28.83 M allocations: 3.707 GiB, 39.39% gc time)
# 8.203 s (28845630 allocations: 3.71 GiB)




# ## import DB and plot to see || data is fine ||
using CSV, DataFrames
df2 = CSV.read("2Bits_DB.csv")
#
# df[1,:]
#
# size(df)
#
# for i = 1:100
#     row = rand(1:size(df)[1])
#     K = df[row,:].K
#     n = df[row,:].n
#     A = df[row,:].A
#     up = df[row,:].up
#     δ =  df[row,:].δ
#     param = [K, n, A, up, δ]
#     println("K:$K n:$n, δ:$δ, A:$A, up:$up")
#     sol, ts = run_prob_2bits(;init_relax = 2000., duration=δ,relax=2000.,signal=A,K=K,n=n, up= up,cycle = 20)
#     plot_2bit(sol, up, param)
# end
##

using VegaLite
using Latexify

eqs2_forlatex = [
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF),
    D(g1)      ~ ξ * hill(p)                   - deg(g1),
    D(g2)      ~ ξ * hill(m1_HKCI)             - deg(g2),
    D(g3)      ~ ξ * hill(g1 + g2)             - deg(g3),
    D(m2_LexA1) ~ ξ * hill(m2_PhlF + g3)       - deg(m2_LexA1),
    D(m2_IcaR ) ~ ξ * hill(m2_LexA1 + g3)      - deg(m2_IcaR),
    D(m2_CI1  ) ~ ξ * hill(m2_LexA1 + m2_PhlF) - deg(m2_CI1),
    D(m2_PsrA ) ~ ξ * hill(m2_IcaR + m2_CI1)   - deg(m2_PsrA),
    D(m2_BM3RI) ~ ξ * hill(m2_PsrA)            - deg(m2_BM3RI),
    D(m2_HKCI ) ~ ξ * hill(m2_BM3RI + m2_PhlF) - deg(m2_HKCI),
    D(m2_PhlF ) ~ ξ * hill(m2_PsrA + m2_HKCI)  - deg(m2_PhlF)]

latexify(eqs2_forlatex)


df2 |> @vlplot(
    :point,
    x=:K,
    y=:n,
    size= "δ:n",
    color=:δ,
    width=400,
    height=400
)
