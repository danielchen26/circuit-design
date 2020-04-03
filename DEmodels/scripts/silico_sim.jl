
#Author: Tianchi Chen
## ------ Import package and functions
using DifferentialEquations, ParameterizedFunctions
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
using Latexify, Random, Base
using CSV, DataFrames, Images
using Parameters, ProgressMeter
using ChemometricsTools
include("functions.jl")
##

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
sol0 = solve(prob0, Tsit5())
plot(sol0, lw = 2, ylims = (0, 22))


#  ==== Induced by signals with t_on time ======
# when constant signal is applied, we expected to see oscillation
p = 20.0;
t_on = 1000.0;
prob1 = ODEProblem(silico, sol0[end], (0.0, t_on), p)
sol1 = solve(prob1, SSRootfind())
plot(sol1, vars = [:m_HKCI, :m_PhlF], lw = 2)







# ======= control for multiple cycles ==============÷
@show Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(δ = 270.0, Δ = 1500.)
prob0 = ODEProblem(silico, u0, tspan, p)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
py = plot(sol, vars = [:m_HKCI, :m_PhlF], lw = 2, xlabel = "time", ylabel = "concentration", title = "Signal Duration: $δ")


tt = zero(ts)
tt[[1,2]] .= 10
tt
scatter!(py,ts,tt)


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







# ======= 1 bit counter =================
# ======= 🔴 TEST:  💚Global optimization ===============
# using BlackBoxOptim, LinearAlgebra
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
function run_prob(;duration,relax)
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(δ = duration, Δ = relax)
    prob0 = ODEProblem(silico, u0, tspan, p)
    sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-15, abstol = 1e-19)
end


# sol = run_prob(duration = δ, relax = 1500.)
# plot(sol, vars = [:m_HKCI, :m_PhlF])


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








## ======= Build multiple counter connectors : 2 Bits counter case 📗 =========

# ==== Define ODEProblem =======
γ = 0.025
ξ = 0.025
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
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
    dg3   = ξ * hill(g1 + g2) - degradation(g3) # g3 sserves as the input for the 2nd bit
    # Bit 2 =============
    dm2_LexA1 = ξ * hill(m2_PhlF + g3) - degradation(m2_LexA1)
    dm2_IcaR  = ξ * hill(m2_LexA1 + g3) - degradation(m2_IcaR)
    dm2_CI1   = ξ * hill(m2_LexA1 + m2_PhlF) - degradation(m2_CI1)
    dm2_PsrA  = ξ * hill(m2_IcaR + m2_CI1) - degradation(m2_PsrA)
    dm2_BM3RI = ξ * hill(m2_PsrA) - degradation(m2_BM3RI)
    dm2_HKCI  = ξ * hill(m2_BM3RI + m2_PhlF) - degradation(m2_HKCI)
    dm2_PhlF  = ξ * hill(m2_PsrA + m2_HKCI) - degradation(m2_PhlF)
end p

##
u0 = Float64[i for i in rand(1:22, 17)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 270, cycle = 10)
prob0 = ODEProblem(silico2, u0, tspan, p)
sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
Ylim = 3*up
C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),size=(600,500))
##
# savefig(Bit2_plt,"./Counter_Plot/2bit_abnormal_case.png")


## function Wraping the above to avoid DomainError
function Bit2_gen(;duration, relE, absE)
    global warn = nothing
    try
        u0 = Float64[i for i in rand(1:22, 17)]
        Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 4000., Δ = 4000., δ = duration, cycle = 10)
        @show tspan
        init_param = Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p
        prob0 = ODEProblem(silico2, u0, tspan, p)
        sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = relE, abstol = absE)
        return sol, init_param, warn
    catch err
        if isa(err, DomainError)
            @show  warn = "DomainError"
        end
    end
    return sol, init_param, warn
end

function Bit2_plot(sol)
    Ylim = 3*up
    C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),size=(600,500))
    display(Bit2_plt)
end
##
sol, init_param, warn = Bit2_gen(duration = 270., relE = 1e-15, absE = 1e-20)
Bit2_plot(sol)
# g6i_av, g6ii_av, g7i_av, g7ii_av = switching(sol,4)
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

## -------- Finding local minmax
C_id = 10 # connector ID in the Diff equation
opt = L_max(sol,C_id, ts[1],ts[2])
py = plot(sol, vars=(0,C_id), plotdensity=10000);
scatter!(py, [opt.minimizer],[-opt.minimum],label="Local Max")
##



##
# -------- First g3 peak detection
C_id = 10 # connector ID in the Diff equation
opt2 = L_max(sol,C_id, ts[3],ts[4])
opt2.minimizer, -opt2.minimum
# --------- Attenuation averaging
function attenuation(sol ;thred = 0.5)
    if -opt2.minimum > thred
        lmax_set = []
        for i = collect(3:4:2*cycle)
            opt_i = L_max(sol,C_id, ts[i],ts[i+1])
            lmax = -opt_i.minimum
            push!(lmax_set,lmax)
        end
        return sum(lmax_set)/length(lmax_set), std(lmax_set)
    end
end

AT_mean, AT_var = attenuation(sol)


function g3_peak_init(sol; ts = ts, thred = 0.5, C_id = 10)
    opt = L_max(sol,C_id, ts[1],ts[2])
    opt2 = L_max(sol,C_id, ts[3],ts[4])
    lmax1 = -opt.minimum;    lmax2 = -opt2.minimum;
    if lmax1 > thred && lmax2 < thred
        g3_peak_i = 1
    elseif lmax2 > thred && lmax1 < thred
        g3_peak_i = 3
    end
end


cycle_id = g3_peak_init(sol; thred = 0.5, C_id = 10)


for ti = 1:2:length(ts)
    @show ti, ti+1
    opt = L_max(sol,C_id, ts[ti],ts[ti + 1])
    lmax1 = -opt.minimum;
    @show lmax1
end



## === Check whether the gate switches or not ============
function switching2(sol, f)
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



























## ------ Import package and functions
using DifferentialEquations, ParameterizedFunctions
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
using Latexify, Random, Base
using CSV, DataFrames, Images
using Parameters, ProgressMeter
using ChemometricsTools
include("functions.jl")

## ==== Build multiple counter connectors : 3 Bits counter case 📗 =========
# ==== Define ODEProblem =======
γ = 0.025
ξ = 0.025
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
degradation(x) = γ * x
up = 1.5;
dn = 0.002;
K = 0.081;
n = 2.81;
silico3 = @ode_def_bare bit3 begin
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
    dg3   = ξ * hill(g1 + g2) - degradation(g3) # g3 sserves as the input for the 2nd bit
    # Bit 2 =============
    dm2_LexA1 = ξ * hill(m2_PhlF + g3) - degradation(m2_LexA1)
    dm2_IcaR  = ξ * hill(m2_LexA1 + g3) - degradation(m2_IcaR)
    dm2_CI1   = ξ * hill(m2_LexA1 + m2_PhlF) - degradation(m2_CI1)
    dm2_PsrA  = ξ * hill(m2_IcaR + m2_CI1) - degradation(m2_PsrA)
    dm2_BM3RI = ξ * hill(m2_PsrA) - degradation(m2_BM3RI)
    dm2_HKCI  = ξ * hill(m2_BM3RI + m2_PhlF) - degradation(m2_HKCI)
    dm2_PhlF  = ξ * hill(m2_PsrA + m2_HKCI) - degradation(m2_PhlF)
    # Connector 2 (Two AND gates) =============
    # 1rst AND gate combines (out1,out2)
    dg21 = ξ * hill(m_HKCI) - degradation(g21)
    dg22  = ξ * hill(m2_HKCI) - degradation(g22)
    dg23   = ξ * hill(g21 + g22) - degradation(g23) # g3 sserves as the input for the 2nd bit
    # 2nd AND gate combines (out g23,p)
    dg24 = ξ * hill(p) - degradation(g24)
    dg25  = ξ * hill(g23) - degradation(g25)
    dg26   = ξ * hill(g24 + g25) - degradation(g26) # g3 sserves as the input for the 2nd bit
    # Bit 3 =============
    dm3_LexA1 = ξ * hill(m3_PhlF + g26) - degradation(m3_LexA1)
    dm3_IcaR  = ξ * hill(m3_LexA1 + g26) - degradation(m3_IcaR)
    dm3_CI1   = ξ * hill(m3_LexA1 + m3_PhlF) - degradation(m3_CI1)
    dm3_PsrA  = ξ * hill(m3_IcaR + m3_CI1) - degradation(m3_PsrA)
    dm3_BM3RI = ξ * hill(m3_PsrA) - degradation(m3_BM3RI)
    dm3_HKCI  = ξ * hill(m3_BM3RI + m3_PhlF) - degradation(m3_HKCI)
    dm3_PhlF  = ξ * hill(m3_PsrA + m3_HKCI) - degradation(m3_PhlF)
end p

##
u0 = Float64[i for i in rand(1:22, 30)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 270, cycle = 20)
prob0 = ODEProblem(silico3, u0, tspan, p)
sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
Ylim = 3*up
# plt_opt = [legend = :topright, legendfontsize = 7, ylims =(0.,Ylim)]
C_plt = plot(sol, vars = [:g3], legend = false, ylims =(0.,Ylim))
C2_plt = plot(sol, vars = [:g26], legend = false, ylims =(0.,Ylim))
B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = false, ylims =(0.,Ylim))
B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = false, ylims =(0.,Ylim))
B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = false, ylims =(0.,Ylim))
Bit2_plt = plot(C_plt,C2_plt,B1_plt,B2_plt,B3_plt,layout = (5,1),size=(600,500))
## optimization for parameters searching

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


function MB_T0(ts_ID, C1, C2; thred =up/2)
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
T0 = MB_T0(ts_ID, C1, C2; thred =up/2)

# function MB_cost(T0, sol; thred =0.5)
#     println("Initial T0:\n", T0)
# end
#
# MB_cost(T0, sol)

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

# Switch(sol,[6,7], Δ0, ts, T0)
# Switch(sol,[6,7], Δ0, ts, T0.+2)
# Switch(sol,[6,7], Δ0, ts, T0.+4)
# Switch(sol,[6,7], Δ0, ts, T0.+6)
#
# Switch(sol,[16,17], Δ0, ts, T0)
# Switch(sol,[16,17], Δ0, ts, T0.+2)
# Switch(sol,[16,17], Δ0, ts, T0.+4)
# Switch(sol,[16,17], Δ0, ts, T0.+6)
#
# Switch(sol,[29,30], Δ0, ts, T0)
# Switch(sol,[29,30], Δ0, ts, T0.+2)
# Switch(sol,[29,30], Δ0, ts, T0.+4)
# Switch(sol,[29,30], Δ0, ts, T0.+6)
# Switch(sol,[29,30], Δ0, ts, T0.+8)
# Switch(sol,[29,30], Δ0, ts, T0.+10)

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

G6,  G7 = Switch_cost(sol, [6,7], Δ0, ts, T0)
G16, G17 = Switch_cost(sol, [16,17], Δ0, ts, T0)
G29, G30 = Switch_cost(sol, [29,30], Δ0, ts, T0)

# --- Make binary number
thred = up/2

b11 = G6 .> thred
b12 = G7 .> thred
@show b11, b12

b21 = G16 .> thred
b22 = G17 .> thred
@show b21, b22

b31 = G29 .> thred
b32 = G30 .> thred
@show b31, b32

# we need to make the below thred gate to be 0, and above the thred gate to be 1
@show b11
@show b21
@show b32



## run ODEProblem for 3 bits counter
# ========  To run a problem 1 time, given input of signal duration
function run_prob_3bits(;duration,relax,signal)
    u0 = Float64[i for i in rand(1:22, 30)]
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = relax, δ = duration, A = signal, cycle = 20)
    prob0 = ODEProblem(silico3, u0, tspan, p)
    sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
end

Ylim = 3*up
δ = 270.
sol = run_prob_3bits(duration = δ, relax = 20000.)
plot(sol, vars = [:m_HKCI, :m_PhlF, :m2_HKCI, :m2_PhlF, :m3_HKCI, :m3_PhlF ], legend = false, layout=(3,2), ylims =(0.,Ylim), title=["Bit1" "Bit1" "Bit2" "Bit2" "Bit3" "Bit3"], title_location=:left)

function cost_bit3(sol, ts)
    # set T0 for multibit counter
    ts_ID, C1, C2  = Carrying_bit(sol, ts)
    T0 = MB_T0(ts_ID, C1, C2; thred =up/2)
    # @show T0
    G6,  G7 = Switch_cost(sol, [6,7], Δ0, ts, T0)
    G16, G17 = Switch_cost(sol, [16,17], Δ0, ts, T0)
    G29, G30 = Switch_cost(sol, [29,30], Δ0, ts, T0)
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
    return costtot
end

costtot = cost_bit3(sol, ts)





# ======== Sampling Parameters ===================
using ProgressMeter
# Varying K, n, δ
df = DataFrame(K = Float64[], n = Float64[], δ =  Float64[])
tt = []
@time @showprogress for K = 0.001: 0.01:0.1, n = 1.5:0.1:3., δ = 250:5:350, A = 5:2:40
    K = K; n = n;
    sol = run_prob_3bits(duration = δ, relax = 20000., signal = A)
    # display(plot(sol, vars = [:m_HKCI, :m_PhlF]))
    # Check switching =============
    costtot == 8 ? push!(df, [K, n, δ]) : nothing

    push!(tt,1.)
    @show sum(tt)
    # if sum(tt) >= 20.
    #     break
    # end
end


# up = 1.5;
# dn = 0.002;
# K = 0.081;
# n = 2.81;

df































## ==== Build multiple counter connectors : 4 Bits counter case 📗 =========
# ==== Define ODEProblem =======
γ = 0.025
ξ = 0.025
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
degradation(x) = γ * x
up = 1.5;
dn = 0.002;
K = 0.081;
n = 2.81;
silico4 = @ode_def_bare bit4 begin
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
    dg3   = ξ * hill(g1 + g2) - degradation(g3) # g3 sserves as the input for the 2nd bit
    # Bit 2 =============
    dm2_LexA1 = ξ * hill(m2_PhlF + g3) - degradation(m2_LexA1)
    dm2_IcaR  = ξ * hill(m2_LexA1 + g3) - degradation(m2_IcaR)
    dm2_CI1   = ξ * hill(m2_LexA1 + m2_PhlF) - degradation(m2_CI1)
    dm2_PsrA  = ξ * hill(m2_IcaR + m2_CI1) - degradation(m2_PsrA)
    dm2_BM3RI = ξ * hill(m2_PsrA) - degradation(m2_BM3RI)
    dm2_HKCI  = ξ * hill(m2_BM3RI + m2_PhlF) - degradation(m2_HKCI)
    dm2_PhlF  = ξ * hill(m2_PsrA + m2_HKCI) - degradation(m2_PhlF)
    # Connector 2 (Two AND gates) =============
    # 1rst AND gate combines (out1,out2)
    dg21 = ξ * hill(m_HKCI) - degradation(g21)
    dg22  = ξ * hill(m2_HKCI) - degradation(g22)
    dg23   = ξ * hill(g21 + g22) - degradation(g23) # g3 sserves as the input for the 2nd bit
    # 2nd AND gate combines (out g23,p)
    dg24 = ξ * hill(p) - degradation(g24)
    dg25  = ξ * hill(g23) - degradation(g25)
    dg26   = ξ * hill(g24 + g25) - degradation(g26) # g3 sserves as the input for the 2nd bit
    # Bit 3 =============
    dm3_LexA1 = ξ * hill(m3_PhlF + g26) - degradation(m3_LexA1)
    dm3_IcaR  = ξ * hill(m3_LexA1 + g26) - degradation(m3_IcaR)
    dm3_CI1   = ξ * hill(m3_LexA1 + m3_PhlF) - degradation(m3_CI1)
    dm3_PsrA  = ξ * hill(m3_IcaR + m3_CI1) - degradation(m3_PsrA)
    dm3_BM3RI = ξ * hill(m3_PsrA) - degradation(m3_BM3RI)
    dm3_HKCI  = ξ * hill(m3_BM3RI + m3_PhlF) - degradation(m3_HKCI)
    dm3_PhlF  = ξ * hill(m3_PsrA + m3_HKCI) - degradation(m3_PhlF)
    # Connector 3 (Two AND gates) =============
    # 1rst AND gate combines (out1,out2)
    dg31 = ξ * hill(m2_HKCI) - degradation(g31)
    dg32  = ξ * hill(m3_HKCI) - degradation(g32)
    dg33   = ξ * hill(g31 + g32) - degradation(g33) # g3 sserves as the input for the 2nd bit
    # 2nd AND gate combines (out g33,p)
    dg34 = ξ * hill(p) - degradation(g34)
    dg35  = ξ * hill(g33) - degradation(g35)
    dg36   = ξ * hill(g34 + g35) - degradation(g36) # g3 sserves as the input for the 2nd bit
    # Bit 4 =============
    dm4_LexA1 = ξ * hill(m4_PhlF + g36) - degradation(m4_LexA1)
    dm4_IcaR  = ξ * hill(m4_LexA1 + g36) - degradation(m4_IcaR)
    dm4_CI1   = ξ * hill(m4_LexA1 + m4_PhlF) - degradation(m4_CI1)
    dm4_PsrA  = ξ * hill(m4_IcaR + m4_CI1) - degradation(m4_PsrA)
    dm4_BM3RI = ξ * hill(m4_PsrA) - degradation(m4_BM3RI)
    dm4_HKCI  = ξ * hill(m4_BM3RI + m4_PhlF) - degradation(m4_HKCI)
    dm4_PhlF  = ξ * hill(m4_PsrA + m4_HKCI) - degradation(m4_PhlF)
end p

##
u0 = Float64[i for i in rand(1:22, 43)]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 2000., Δ = 2000., δ = 275, cycle = 20)
prob0 = ODEProblem(silico4, u0, tspan, p)
sol = solve(prob0, SSRootfind(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
Ylim = 3*up
C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
C2_plt = plot(sol, vars = [:g26], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
C3_plt = plot(sol, vars = [:g36], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B1_plt = plot(sol, vars = [:m_HKCI, :m_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B3_plt = plot(sol, vars = [:m3_HKCI, :m3_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
B4_plt = plot(sol, vars = [:m4_HKCI, :m4_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
Bit2_plt = plot(C_plt,C2_plt,C3_plt,B1_plt,B2_plt,B3_plt,B4_plt,layout = (7,1),size=(600,500))
##

using Blink
w = Window()
body!(w, Bit2_plt)




##

time, signal = signal_gen(10,2000.,2000.,270.,20.)
##
