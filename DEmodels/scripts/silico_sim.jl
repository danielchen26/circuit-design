#Author: Tianchi Chen
using DifferentialEquations, ParameterizedFunctions
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
using Latexify, Random, Base
using CSV, DataFrames, Images
using Parameters, ProgressMeter
include("functions.jl")

using BlackBoxOptim, LinearAlgebra

γ = 0.025
ξ = 0.025
hill(x) = dn + (up - dn) * K^n / (K^n + x^n)
degradation(x) = γ * x

up = 1.5; dn = 0.002; K = 0.081; n = 2.81
silico = @ode_def_bare counter begin
    dm_LexA1 = ξ * hill(m_PhlF + p) -  degradation(m_LexA1)
    dm_IcaR  = ξ * hill(m_LexA1 + p) - degradation(m_IcaR)
    dm_CI1   = ξ * hill(m_LexA1 + m_PhlF) - degradation(m_CI1)
    dm_PsrA  = ξ * hill(m_IcaR + m_CI1) - degradation(m_PsrA)
    dm_BM3RI = ξ * hill(m_PsrA) - degradation(m_BM3RI)
    dm_HKCI  = ξ * hill(m_BM3RI + m_PhlF) - degradation(m_HKCI)
    dm_PhlF  = ξ * hill(m_PsrA + m_HKCI) - degradation(m_PhlF)
end p


# ==== Randomized Initials equlibrations =====
u0 = Float64[i for i in rand(1:22, 7)]
p = 0.0
tspan = (0.0, 3000.0)
prob0 = ODEProblem(silico, u0, tspan, p)
sol0 = solve(prob0, SSRootfind())
plot(sol0, lw = 2, ylims = (0, 22))


#  ==== Induced by signals with t_on time ======
# when constant signal is applied, we expected to see oscillation
p = 20.0;
t_on = 1000.0;
prob1 = ODEProblem(silico, sol0[end], (0.0, t_on), p)
sol1 = solve(prob1, SSRootfind())
plot(sol1, vars = [:m_HKCI, :m_PhlF], lw = 2)







# ======= control for multiple cycles ==============÷
@show Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(δ = 270.0)
prob0 = ODEProblem(silico, u0, tspan, p)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
py = plot(sol, vars = [:m_HKCI, :m_PhlF], lw = 2, xlabel = "time", ylabel = "concentration", title = "Signal Duration: $δ")


tt = zero(ts)
tt[[1,2]] .= 10
tt
scatter(ts,tt)



# ===== Animation of control problems =======
sig_rg = 250:10:330
anim = Anim_gen_controlproblem(cycle, Δ0, Δ, δ, A, sig_rg)
gif(anim, "/tmp/anim.gif", fps = 5)






# Find the lb and ub of the P inpluse and the min time to wait for the next signal
# Find the local minimum index
#  ========== Mark the local minimum ==========

# ====== Local minimum and maximum with t_on ============
t_on = 1000.0
plt_min, mark6_min, mark7_min = localminmax_sig(sol1, t_on, "lmin")
plt_min
plt_max, mark6_max, mark7_max = localminmax_sig(sol1, t_on, "lmax")
plt_max







#  Test the time starting from ti =1 to ti = first local min
# t_wait = []
ul_range = []
t_range = []
for ti = 1:mark7_min[2]
    u0_t11 = sol1[:, ti]
    println("Signal Duration : ", sol1.t[ti])
    p = 0
    prob11 = ODEProblem(silico, u0_t11, (0.0, 3000.0), p)
    sol11 = solve(prob11, Tsit5())
    display(plot(sol11, vars = [(0, 6), (0, 7)]))

    # How much time to wait until the system become statble again
    locs = findlocalmaxima(sol11[7, :])#[1]
    marklist = [i[1] for i in locs]
    marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
    display(scatter!(sol11[marklist], vars = [:m_PhlF], xlims = (1:t_on), marker = marker_spec))
    # display(plot(sol11[6,:]))
    # stable_t_ind = locs[1]
    if sol11[7, end] > sol11[6, end]
        push!(ul_range, ti)
        push!(t_range, sol1.t[ti])
        # push!(t_wait, sol11.t[stable_t_ind])
    end
end

@show t_range
@show ul_range


ti = 45
u0_t11 = sol1[:, ti]
println("Signal Duration : ", sol1.t[ti])
p = 0
prob11 = ODEProblem(silico, u0_t11, (0.0, 3000.0), p)
sol11 = solve(prob11, Tsit5())
display(plot(sol11, vars = [(0, 6), (0, 7)]))

locs = findlocalmaxima(sol11[7, :])#[1]
marklist = [i[1] for i in locs]
marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
display(scatter!(sol11[marklist], vars = [:m_PhlF], xlims = (1:t_on), marker = marker_spec))





#  give the min time to wait(lower bound) for the next signal within the switching range
maximum(t_wait)
# switching index range
ul_range
# switching lower bound time
t_lb = t_range[1]
# switching upper bound time
t_ub = t_range[end]








# ======= 🔴 TEST:  💚Global optimization ===============
using BlackBoxOptim, LinearAlgebra
# plotly()
# ==== Define ODEProblem =======
γ = 0.025
ξ = 0.025
hill(x) = dn + (up - dn) * K^n / (K^n + x^n)
degradation(x) = γ * x
up = 1.5;
dn = 0.002;
K = 0.081;
n = 2.81;
silico = @ode_def_bare counter begin
    dm_LexA1 = ξ * hill(m_PhlF + p) - degradation(m_LexA1)
    dm_IcaR  = ξ * hill(m_LexA1 + p) - degradation(m_IcaR)
    dm_CI1   = ξ * hill(m_LexA1 + m_PhlF) - degradation(m_CI1)
    dm_PsrA  = ξ * hill(m_IcaR + m_CI1) - degradation(m_PsrA)
    dm_BM3RI = ξ * hill(m_PsrA) - degradation(m_BM3RI)
    dm_HKCI  = ξ * hill(m_BM3RI + m_PhlF) - degradation(m_HKCI)
    dm_PhlF  = ξ * hill(m_PsrA + m_HKCI) - degradation(m_PhlF)
end p



u0 = Float64[i for i in rand(1:22, 7)]
@show Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(δ = 270, cycle = 4)
prob0 = ODEProblem(silico, u0, tspan, p)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
plot(sol, vars = [:m_HKCI, :m_PhlF])



function cost(pa)
    # try to optimize δ: signal duration, K: dissociation rate, n: hill para
    up = pa[1]
    K = pa[2]
    n = pa[3]
    @show Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control()# δ = param
    prob0 = ODEProblem(silico, u0, tspan, p)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
    norm(sol[7, end] - up)
end

bound = [(1.5, 5.5), (1.5, 5.5), (1.5, 5.5)]
result = bboptimize(cost; SearchRange = bound, NumDimensions = 3)





function run_prob(;duration)
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(δ = duration)
    prob0 = ODEProblem(silico, u0, tspan, p)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
end




# ======== Sampling Parameters ===================
using ProgressMeter
# Varying K, n, δ
df = DataFrame(K = Float64[], n = Float64[], δ =  Float64[])
tt = []
@time @showprogress for K = 0.001: 0.01:0.1, n = 1.5:0.1:3., δ = 250:5:350
    K = K; n = n;
    sol = run_prob(duration = δ)
    # display(plot(sol, vars = [:m_HKCI, :m_PhlF]))
    # Check switching =============
    g6i_av, g6ii_av, g7i_av, g7ii_av = switching(sol,4)

    if g6i_av > g7i_av
        dn < g6ii_av < (up+dn)/2  &&  (up+dn)/2 < g7ii_av < up ? push!(df, [K, n, δ]) : nothing
    elseif g6i_av < g7i_av
        dn < g7ii_av < (up+dn)/2  &&   (up+dn)/2 < g6ii_av < up ? push!(df, [K, n, δ]) : nothing
    end

    push!(tt,1.)
    @show sum(tt)
    # if sum(tt) >= 200.
    #     break
    # end
end
df
CSV.write("counter_db.csv", df)




# ======= Can make an animation =========
# @time for i = 1: size(df,1)
#     K = df.K[i]; n= df.n[i]; δ =  df.δ[i]
#     sol = run_prob(duration = δ)
#     # display(plot(sol, vars = [:m_HKCI, :m_PhlF]))
# end








# ======= Build multiple counter connectors : 2 Bits counter case 📗 =========

# ==== Define ODEProblem =======
γ = 0.025
ξ = 0.025
hill(x) = dn + (up - dn) * K^n / (K^n + x^n)
degradation(x) = γ * x
up = 1.5;
dn = 0.002;
K = 0.081;
n = 2.81;
silico2 = @ode_def_bare bit2 begin
    # Bit 1 =================
    dm_LexA1 = ξ * hill(m_PhlF + p) - degradation(m_LexA1)
    dm_IcaR  = ξ * hill(m_LexA1 + p) - degradation(m_IcaR)
    dm_CI1   = ξ * hill(m_LexA1 + m_PhlF) - degradation(m_CI1)
    dm_PsrA  = ξ * hill(m_IcaR + m_CI1) - degradation(m_PsrA)
    dm_BM3RI = ξ * hill(m_PsrA) - degradation(m_BM3RI)
    dm_HKCI  = ξ * hill(m_BM3RI + m_PhlF) - degradation(m_HKCI)
    dm_PhlF  = ξ * hill(m_PsrA + m_HKCI) - degradation(m_PhlF)
    # Connector 1 =============
    dg1 = ξ * hill(p) - degradation(g1)
    dg2  = ξ * hill(m_HKCI) - degradation(g2)
    dg3   = ξ * hill(g1 + g2) - degradation(g3)
    # Bit 2 =============
    dm2_LexA1 = ξ * hill(m2_PhlF + g3) - degradation(m2_LexA1)
    dm2_IcaR  = ξ * hill(m2_LexA1 + g3) - degradation(m2_IcaR)
    dm2_CI1   = ξ * hill(m2_LexA1 + m2_PhlF) - degradation(m2_CI1)
    dm2_PsrA  = ξ * hill(m2_IcaR + m2_CI1) - degradation(m2_PsrA)
    dm2_BM3RI = ξ * hill(m2_PsrA) - degradation(m2_BM3RI)
    dm2_HKCI  = ξ * hill(m2_BM3RI + m2_PhlF) - degradation(m2_HKCI)
    dm2_PhlF  = ξ * hill(m2_PsrA + m2_HKCI) - degradation(m2_PhlF)
end p


u0 = Float64[i for i in rand(1:22, 17)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 270, cycle = 10)
prob0 = ODEProblem(silico2, u0, tspan, p)
sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)

Ylim = 3*up
C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),size=(600,500))

# savefig(Bit2_plt,"./Counter_Plot/2bit_abnormal_case.png")





u0 = Float64[i for i in rand(1:22, 17)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 270, cycle = 10)
prob0 = ODEProblem(silico2, u0, tspan, p)
sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)

C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,5.))
B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,5.))
B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,5.))
Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),size=(600,500))



function Bit2_gen(;duration, relE, absE)
    global warn = nothing
    try
        u0 = Float64[i for i in rand(1:22, 17)]
        Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = duration, cycle = 10)
        prob0 = ODEProblem(silico2, u0, tspan, p)
        sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = relE, abstol = absE)

        C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,5.))
        B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,5.))
        B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,5.))
        Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),size=(600,500))
        display(Bit2_plt)
    catch err
        if isa(err, DomainError)
            @show  warn = "DomainError"
        end
    end
    return sol, warn
end

function Bit2_gen_F(;duration, relE, absE) # if domain error exist, try other initial condition
    sol, warn = Bit2_gen(duration = duration, relE = relE, absE = absE)
    while warn == "DomainError"
        sol, warn = Bit2_gen(duration = 270., relE = 1e-15, absE = 1e-20)
    end
    return sol, warn
end

sol, warn = Bit2_gen_F(duration = 270., relE = 1e-15, absE = 1e-20)
g6i_av, g6ii_av, g7i_av, g7ii_av = switching(sol,4)


using ChemometricsTools
cc = DataFrame(a=sol[10,:])
CSV.write("cc.csv",cc)
cc = CSV.read("ccc.csv")
plot(cc.a)
locs = findpeaks(cc.a, m=200)#findlocalmaxima(cc.a)
marklist = [i[1] for i in locs]
marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
# scatter!(sol[marklist], vars = [:g3], xlims = (1:2000.0), marker = marker_spec) # xlims = (1:t_on),
scatter!(marklist, cc.a[marklist],marker = marker_spec) # xlims = (1:t_on),

using ChemometricsTools
cc
plot(cc.a)
findpeaks(cc.a, m=10)

Y =collect(cc.a)
filtered_Y = Y
filtered_Y[ Y .< 1] .= 0.0


locs = findpeaks(filtered_Y; m = 3)
marklist = [i[1] for i in locs]
marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
# scatter!(sol[marklist], vars = [:g3], xlims = (1:2000.0), marker = marker_spec) # xlims = (1:t_on),
scatter!(marklist, cc.a[marklist],marker = marker_spec) # xlims = (1:t_on),



using OrdinaryDiffEq
initial = [0.01, 0.01, 0.01, 0.01]
tspan = (0.,100.)

#Define the problem
function double_pendulum_hamiltonian(udot,u,p,t)
    α  = u[1]
    lα = u[2]
    β  = u[3]
    lβ = u[4]
    udot .=
    [2(lα-(1+cos(β))lβ)/(3-cos(2β)),
    -2sin(α) - sin(α+β),
    2(-(1+cos(β))lα + (3+2cos(β))lβ)/(3-cos(2β)),
    -sin(α+β) - 2sin(β)*(((lα-lβ)lβ)/(3-cos(2β))) + 2sin(2β)*((lα^2 - 2(1+cos(β))lα*lβ + (3+2cos(β))lβ^2)/(3-cos(2β))^2)]
end

#Pass to solvers
poincare = ODEProblem(double_pendulum_hamiltonian, initial, tspan)
sol = solve(poincare, Tsit5())
f = (t) -> -sol(first(t),idxs=4)

lower = [0.02]
upper = [20.03]
initial_x = [15.02]
inner_optimizer = LBFGS()
opt2 = optimize(f, lower, upper, initial_x, Fminbox(inner_optimizer))

plot(sol, vars=(0,4), plotdensity=10000)
scatter!([opt2.minimizer],[-opt2.minimum],label="Local Max")






# === Check whether the gate switches or not ============
function switching2(sol, cycles)
    # println("End of relaxation time: ", Δ0)
    G6 = sol(Δ0)[6];  G7 = sol(Δ0)[7]
    G62 = sol(Δ0)[16];  G72 = sol(Δ0)[17]
    @show G6, G7
    @show G62, G72
    g6ii_tot = g7ii_tot = 0
    g6i_tot = g7i_tot = 0
    if G6 > G7

        for i = 1: 2: cycles
            tii = Δ0 + i * Δ
            ti = Δ0 + (i - 1) * Δ
            # @show ti, tii
            g6ii = sol(tii)[6];    g6i = sol(ti)[6]
            g7ii = sol(tii)[7];    g7i = sol(ti)[7]
            # @show g6i, g6ii
            # @show g7i, g7ii

            g6ii_tot += g6ii; g7ii_tot += g7ii
            g6i_tot += g6i  ; g7i_tot += g7i
        end
    elseif G6 < G7
        for i = 1: 2: cycles
            tii = Δ0 + i * Δ
            ti = Δ0 + (i - 1) * Δ
            # @show ti, tii
            g6ii = sol(tii)[6];    g6i = sol(ti)[6]
            g7ii = sol(tii)[7];    g7i = sol(ti)[7]
            # @show g6i, g6ii
            # @show g7i, g7ii

            g6ii_tot += g6ii; g7ii_tot += g7ii
            g6i_tot += g6i  ; g7i_tot += g7i
        end
    end
    cnt = length(collect(1: 2: cycles))
    g6i_av = g6i_tot/cnt; g6ii_av = g6ii_tot/cnt
    g7i_av = g7i_tot/cnt; g7ii_av = g7ii_tot/cnt
    # @show g6i_av, g6ii_av
    # @show g7i_av, g7ii_av
    return g6i_av, g6ii_av, g7i_av, g7ii_av
end




using ProgressMeter
# Varying K, n, δ
df = DataFrame(K = Float64[], n = Float64[], δ =  Float64[])
tt = []
@time @showprogress for K = 0.001: 0.01:0.1, n = 1.5:0.1:3., δ = 250:5:350
    K = K; n = n;
    sol, warn = Bit2_gen_F(duration = δ,relE = 1e-15, absE = 1e-20)

    # Check switching =============
    g6i_av, g6ii_av, g7i_av, g7ii_av = switching(sol,4)

    if g6i_av > g7i_av
        dn < g6ii_av < (up+dn)/2  &&  (up+dn)/2 < g7ii_av < up ? push!(df, [K, n, δ]) : nothing
    elseif g6i_av < g7i_av
        dn < g7ii_av < (up+dn)/2  &&   (up+dn)/2 < g6ii_av < up ? push!(df, [K, n, δ]) : nothing
    end

    push!(tt,1.)
    @show sum(tt)
    if sum(tt) >= 30.
        break
    end
end


df






# ========= Cost function for 2 bits ===================
function cost2
end


function run_prob(;prob, duration)
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(δ = duration)
    prob0 = ODEProblem(prob, u0, tspan, p)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
end

# ======== Sampling Parameters for 2 bits counter ===================

# Varying K, n, δ
df = DataFrame(K = Float64[], n = Float64[], δ =  Float64[])
tt = []
@time @showprogress for K = 0.001: 0.01:0.1, n = 1.5:0.1:3., δ = 250:5:350
    K = K; n = n;
    sol = run_prob(duration = δ)
    # display(plot(sol, vars = [:m_HKCI, :m_PhlF]))
    # Check switching =============
    g6i_av, g6ii_av, g7i_av, g7ii_av = switching(sol,4)

    if g6i_av > g7i_av
        dn < g6ii_av < (up+dn)/2  &&  (up+dn)/2 < g7ii_av < up ? push!(df, [K, n, δ]) : nothing
    elseif g6i_av < g7i_av
        dn < g7ii_av < (up+dn)/2  &&   (up+dn)/2 < g6ii_av < up ? push!(df, [K, n, δ]) : nothing
    end

    push!(tt,1.)
    @show sum(tt)
    # if sum(tt) >= 200.
    #     break
    # end
end
