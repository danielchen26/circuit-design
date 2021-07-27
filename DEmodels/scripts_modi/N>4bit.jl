# -*- coding: utf-8 -*-
## Import Functions
using ModelingToolkit, DifferentialEquations, Plots;pyplot()
using Optim
using Printf, DataFrames, CSV
using LaTeXStrings
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




## ==================================================================================================================
## ==================================================================================================================
## ==================================================================================================================

## 4bit counter
## ==== Build multiple counter connectors : 3 Bits counter case 📗 =========
### Define a differential equation system
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
function deg(x)
    γ * x
end
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t) g1(t) g2(t) g3(t)
@variables m2_LexA1(t) m2_IcaR(t) m2_CI1(t) m2_PsrA(t) m2_BM3RI(t) m2_HKCI(t) m2_PhlF(t) g21(t) g22(t) g23(t) g24(t) g25(t) g26(t)
@variables m3_LexA1(t) m3_IcaR(t) m3_CI1(t) m3_PsrA(t) m3_BM3RI(t) m3_HKCI(t) m3_PhlF(t) g31(t) g32(t) g33(t) g34(t) g35(t) g36(t)
@variables m4_LexA1(t) m4_IcaR(t) m4_CI1(t) m4_PsrA(t) m4_BM3RI(t) m4_HKCI(t) m4_PhlF(t) g41(t) g42(t) g43(t) g44(t) g45(t) g46(t)
@variables m5_LexA1(t) m5_IcaR(t) m5_CI1(t) m5_PsrA(t) m5_BM3RI(t) m5_HKCI(t) m5_PhlF(t) g51(t) g52(t) g53(t) g54(t) g55(t) g56(t)
@variables m6_LexA1(t) m6_IcaR(t) m6_CI1(t) m6_PsrA(t) m6_BM3RI(t) m6_HKCI(t) m6_PhlF(t) g61(t) g62(t) g63(t) g64(t) g65(t) g66(t)
@variables m7_LexA1(t) m7_IcaR(t) m7_CI1(t) m7_PsrA(t) m7_BM3RI(t) m7_HKCI(t) m7_PhlF(t) g71(t) g72(t) g73(t) g74(t) g75(t) g76(t)
@variables m8_LexA1(t) m8_IcaR(t) m8_CI1(t) m8_PsrA(t) m8_BM3RI(t) m8_HKCI(t) m8_PhlF(t)
D = Differential(t)

eqs = [
    # 🔵 Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF),

    # 📗 Connector 1 ============= input (p, B1)
    D(g1)      ~ ξ * hill(p)                   - deg(g1),
    D(g2)      ~ ξ * hill(m1_HKCI)             - deg(g2),
    D(g3)      ~ ξ * hill(g1 + g2)             - deg(g3), # g3 serves as the input for the 2nd bit

    # 🔵 Bit 2 ============= input (C1)
    D(m2_LexA1) ~ ξ * hill(m2_PhlF + g3)       - deg(m2_LexA1),
    D(m2_IcaR ) ~ ξ * hill(m2_LexA1 + g3)      - deg(m2_IcaR),
    D(m2_CI1  ) ~ ξ * hill(m2_LexA1 + m2_PhlF) - deg(m2_CI1),
    D(m2_PsrA ) ~ ξ * hill(m2_IcaR + m2_CI1)   - deg(m2_PsrA),
    D(m2_BM3RI) ~ ξ * hill(m2_PsrA)            - deg(m2_BM3RI),
    D(m2_HKCI ) ~ ξ * hill(m2_BM3RI + m2_PhlF) - deg(m2_HKCI),
    D(m2_PhlF ) ~ ξ * hill(m2_PsrA + m2_HKCI)  - deg(m2_PhlF),

    # 📗 Connector 2 (Two AND gates) ============= input (p, C1, B2)
    # 1rst AND gate combines (out1,out2)
    D(g21)     ~ ξ * hill(m1_HKCI)              - deg(g21), #
    D(g22)     ~ ξ * hill(m2_HKCI)             - deg(g22),
    D(g23)     ~ ξ * hill(g21 + g22)           - deg(g23), # g3 serves as the input for the 2nd bit
    # 2nd AND gate combines (out g23,p)
    D(g24)     ~ ξ * hill(p)                   - deg(g24),
    D(g25)     ~ ξ * hill(g23)                 - deg(g25),
    D(g26)     ~ ξ * hill(g24 + g25)           - deg(g26), # g3 serves as the input for the 2nd bit

    # 🔵 Bit 3 ============= input(C2)
    D(m3_LexA1) ~ ξ * hill(m3_PhlF + g26)       - deg(m3_LexA1),
    D(m3_IcaR ) ~ ξ * hill(m3_LexA1 + g26)      - deg(m3_IcaR),
    D(m3_CI1  ) ~ ξ * hill(m3_LexA1 + m3_PhlF)  - deg(m3_CI1),
    D(m3_PsrA ) ~ ξ * hill(m3_IcaR + m3_CI1)    - deg(m3_PsrA),
    D(m3_BM3RI) ~ ξ * hill(m3_PsrA)             - deg(m3_BM3RI),
    D(m3_HKCI ) ~ ξ * hill(m3_BM3RI + m3_PhlF)  - deg(m3_HKCI),
    D(m3_PhlF ) ~ ξ * hill(m3_PsrA + m3_HKCI)   - deg(m3_PhlF),

    # 📗 Connector 3 (Two AND gates) ============= input( p, C2, B3 )
    # 1rst AND gate combines (out1,out2)
    D(g31)     ~ ξ * hill(g23)                 - deg(g31), # 🔴 m2_HKCI changed to g26 (feed forward the output of C2 to C3)
    D(g32)     ~ ξ * hill(m3_HKCI)             - deg(g32),
    D(g33)     ~ ξ * hill(g31 + g32)           - deg(g33), # g3 serves as the input for the 2nd bit
    # 2nd AND gate combines (out g23,p)
    D(g34)     ~ ξ * hill(p)                   - deg(g34),
    D(g35)     ~ ξ * hill(g33)                 - deg(g35),
    D(g36)     ~ ξ * hill(g34 + g35)           - deg(g36), # g3 serves as the input for the 2nd bit

    # 🔵 Bit 4 ============= ininput(C3）
    D(m4_LexA1) ~ ξ * hill(m4_PhlF + g36)       - deg(m4_LexA1),
    D(m4_IcaR ) ~ ξ * hill(m4_LexA1 + g36)      - deg(m4_IcaR),
    D(m4_CI1  ) ~ ξ * hill(m4_LexA1 + m4_PhlF)  - deg(m4_CI1),
    D(m4_PsrA ) ~ ξ * hill(m4_IcaR + m4_CI1)    - deg(m4_PsrA),
    D(m4_BM3RI) ~ ξ * hill(m4_PsrA)             - deg(m4_BM3RI),
    D(m4_HKCI ) ~ ξ * hill(m4_BM3RI + m4_PhlF)  - deg(m4_HKCI),
    D(m4_PhlF ) ~ ξ * hill(m4_PsrA + m4_HKCI)   - deg(m4_PhlF),

        # 📗 Connector 4 (Two AND gates) ============= input( p, C3, B4 )
    # 1rst AND gate combines (out1,out2)
    D(g41)     ~ ξ * hill(g33)                 - deg(g41), # 🔴 m2_HKCI changed to g26 (feed forward the output of C2 to C3)
    D(g42)     ~ ξ * hill(m4_HKCI)             - deg(g42),
    D(g43)     ~ ξ * hill(g41 + g42)           - deg(g43), # g4 serves as the input for the 2nd bit
    # 2nd AND gate combines (out g23,p)
    D(g44)     ~ ξ * hill(p)                   - deg(g44),
    D(g45)     ~ ξ * hill(g43)                 - deg(g45),
    D(g46)     ~ ξ * hill(g44 + g45)           - deg(g46), # g3 serves as the input for the 2nd bit

    # 🔵 Bit 5 ============= ininput(C4）
    D(m5_LexA1) ~ ξ * hill(m5_PhlF + g46)       - deg(m5_LexA1),
    D(m5_IcaR ) ~ ξ * hill(m5_LexA1 + g46)      - deg(m5_IcaR),
    D(m5_CI1  ) ~ ξ * hill(m5_LexA1 + m5_PhlF)  - deg(m5_CI1),
    D(m5_PsrA ) ~ ξ * hill(m5_IcaR + m5_CI1)    - deg(m5_PsrA),
    D(m5_BM3RI) ~ ξ * hill(m5_PsrA)             - deg(m5_BM3RI),
    D(m5_HKCI ) ~ ξ * hill(m5_BM3RI + m5_PhlF)  - deg(m5_HKCI),
    D(m5_PhlF ) ~ ξ * hill(m5_PsrA + m5_HKCI)   - deg(m5_PhlF),

    # 📗 Connector 5 (Two AND gates) ============= input( p, C4, B5 )
    # 1rst AND gate combines (out1,out2)
    D(g51)     ~ ξ * hill(g43)                 - deg(g51), # 🔴
    D(g52)     ~ ξ * hill(m5_HKCI)             - deg(g52),
    D(g53)     ~ ξ * hill(g51 + g52)           - deg(g53), # gN3 serves as the input for the next connector
    # 2nd AND gate combines (out g23,p)
    D(g54)     ~ ξ * hill(p)                   - deg(g54),
    D(g55)     ~ ξ * hill(g53)                 - deg(g55),
    D(g56)     ~ ξ * hill(g54 + g55)           - deg(g56), # gN6 serves as the input for the next bit

    # 🔵 Bit 6 ============= ininput(C5）
    D(m6_LexA1) ~ ξ * hill(m6_PhlF + g56)       - deg(m6_LexA1),
    D(m6_IcaR ) ~ ξ * hill(m6_LexA1 + g56)      - deg(m6_IcaR),
    D(m6_CI1  ) ~ ξ * hill(m6_LexA1 + m6_PhlF)  - deg(m6_CI1),
    D(m6_PsrA ) ~ ξ * hill(m6_IcaR + m6_CI1)    - deg(m6_PsrA),
    D(m6_BM3RI) ~ ξ * hill(m6_PsrA)             - deg(m6_BM3RI),
    D(m6_HKCI ) ~ ξ * hill(m6_BM3RI + m6_PhlF)  - deg(m6_HKCI),
    D(m6_PhlF ) ~ ξ * hill(m6_PsrA + m6_HKCI)   - deg(m6_PhlF),

    # 📗 Connector 6 (Two AND gates) ============= input( p, C5, B6 )
    # 1rst AND gate combines (out1,out2)
    D(g61)     ~ ξ * hill(g53)                 - deg(g61), # 🔴
    D(g62)     ~ ξ * hill(m6_HKCI)             - deg(g62),
    D(g63)     ~ ξ * hill(g61 + g62)           - deg(g63), # gN3 serves as the input for the next connector
    # 2nd AND gate combines (out g23,p)
    D(g64)     ~ ξ * hill(p)                   - deg(g64),
    D(g65)     ~ ξ * hill(g63)                 - deg(g65),
    D(g66)     ~ ξ * hill(g64 + g65)           - deg(g66), # gN6 serves as the input for the next bit

    # 🔵 Bit 7 ============= ininput(C6）
    D(m7_LexA1) ~ ξ * hill(m7_PhlF + g66)       - deg(m7_LexA1),
    D(m7_IcaR ) ~ ξ * hill(m7_LexA1 + g66)      - deg(m7_IcaR),
    D(m7_CI1  ) ~ ξ * hill(m7_LexA1 + m7_PhlF)  - deg(m7_CI1),
    D(m7_PsrA ) ~ ξ * hill(m7_IcaR + m7_CI1)    - deg(m7_PsrA),
    D(m7_BM3RI) ~ ξ * hill(m7_PsrA)             - deg(m7_BM3RI),
    D(m7_HKCI ) ~ ξ * hill(m7_BM3RI + m7_PhlF)  - deg(m7_HKCI),
    D(m7_PhlF ) ~ ξ * hill(m7_PsrA + m7_HKCI)   - deg(m7_PhlF),

    # 📗 Connector 7 (Two AND gates) ============= input( p, C6, B7 )
    # 1rst AND gate combines (out1,out2)
    D(g71)     ~ ξ * hill(g63)                 - deg(g71), # 🔴
    D(g72)     ~ ξ * hill(m7_HKCI)             - deg(g72),
    D(g73)     ~ ξ * hill(g71 + g72)           - deg(g73), # gN3 serves as the input for the next connector
    # 2nd AND gate combines (out g23,p)
    D(g74)     ~ ξ * hill(p)                   - deg(g74),
    D(g75)     ~ ξ * hill(g73)                 - deg(g75),
    D(g76)     ~ ξ * hill(g74 + g75)           - deg(g76), # gN6 serves as the input for the next bit

    # 🔵 Bit 8 ============= ininput(C7）
    D(m8_LexA1) ~ ξ * hill(m8_PhlF + g76)       - deg(m8_LexA1),
    D(m8_IcaR ) ~ ξ * hill(m8_LexA1 + g76)      - deg(m8_IcaR),
    D(m8_CI1  ) ~ ξ * hill(m8_LexA1 + m8_PhlF)  - deg(m8_CI1),
    D(m8_PsrA ) ~ ξ * hill(m8_IcaR + m8_CI1)    - deg(m8_PsrA),
    D(m8_BM3RI) ~ ξ * hill(m8_PsrA)             - deg(m8_BM3RI),
    D(m8_HKCI ) ~ ξ * hill(m8_BM3RI + m8_PhlF)  - deg(m8_HKCI),
    D(m8_PhlF ) ~ ξ * hill(m8_PsrA + m8_HKCI)   - deg(m8_PhlF)
    ]
de = ODESystem(eqs, t,
                [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF, g1, g2, g3,
                 m2_LexA1, m2_IcaR, m2_CI1, m2_PsrA,m2_BM3RI, m2_HKCI, m2_PhlF, g21, g22, g23, g24, g25, g26,
                 m3_LexA1, m3_IcaR, m3_CI1, m3_PsrA, m3_BM3RI, m3_HKCI, m3_PhlF, g31, g32, g33, g34, g35, g36,
                 m4_LexA1, m4_IcaR, m4_CI1, m4_PsrA, m4_BM3RI, m4_HKCI, m4_PhlF, g41, g42, g43, g44, g45, g46,
                 m5_LexA1, m5_IcaR, m5_CI1, m5_PsrA, m5_BM3RI, m5_HKCI, m5_PhlF, g51, g52, g53, g54, g55, g56,
                 m6_LexA1, m6_IcaR, m6_CI1, m6_PsrA, m6_BM3RI, m6_HKCI, m6_PhlF, g61, g62, g63, g64, g65, g66,
                 m7_LexA1, m7_IcaR, m7_CI1, m7_PsrA, m7_BM3RI, m7_HKCI, m7_PhlF, g71, g72, g73, g74, g75, g76,
                 m8_LexA1, m8_IcaR, m8_CI1, m8_PsrA, m8_BM3RI, m8_HKCI, m8_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f = ODEFunction(de)

de
## Solve DifferentialEquations
u0 =  Float64[i for i in rand(1:22, length(eqs))]
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 270, cycle = 240)
param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
prob0 = ODEProblem(ode_f, u0, tspan, param)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)

## Wrap the above problem to a fucntion
function run_prob_5bits(;init_relax,duration,relax,signal,K,n,up,cycle)
    u0 =  Float64[i for i in rand(1:22, length(eqs))]
    Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = init_relax, Δ = relax, δ = duration, A = signal, cycle = cycle)
    param = [up,0.002,K,n,0.025,0.025,p]
    prob0 = ODEProblem(ode_f, u0, (0,init_relax), param)
    sol0 = solve(prob0, Tsit5())
    prob = ODEProblem(ode_f, sol0[end], tspan, param)
    sol = solve(prob, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)
    return sol, ts
end

## ploting
function plot_5bit(sol, ts, up; opt = 'C')
    Ylim = 2*up
    if opt == 'C'
        C_plt = plot(sol, vars = [g3], label = "Connector type 1",
                    legendtitle = "Connector module",
                    # ylims =(0.,Ylim)
                    )
        C2_plt = plot(sol, vars = [g26], label = "Connector type 2",
                    legendtitle = "Connector module",
                    # ylims =(0.,Ylim)
                    )
        C3_plt = plot(sol, vars = [g36], label = "Connector type 2",
                    legendtitle = "Connector module",
                    # ylims =(0.,Ylim)
                    )
        C4_plt = plot(sol, vars = [g46], label = "Connector type 2",
                    legendtitle = "Connector module",
                    # ylims =(0.,Ylim)
                    )
        Bit5_plt_C = plot(
                        C_plt, C2_plt, C3_plt, C4_plt,
                        # B1_plt,C_plt,B2_plt,C2_plt,B3_plt,C3_plt,B4_plt,
                        layout = (4,1),
                        xtickfontsize=11,ytickfontsize=11, size=(1000,800),
                        legend = :topright,
                        # legend = :outertopright,
                        legendfont = font("Times new roman", 10),
                        xlabel = "Time", ylabel = "Concentration",dpi = 400,
                        # title = ["Carrying bit 1" "Carrying bit 2" "B1" "B2" "B3"]
                        # title = ["B1" "Carrying bit 1" "B2" "Carrying bit 2" "B3"]
                        )
        # display(Bit3_plt)
        return Bit5_plt_C

    elseif opt == 'B'
        B1_plt = plot(sol, vars = [m1_HKCI, m1_PhlF],
                    label = ["Q1: output" L"$\bar{Q}$1: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The first bit counter",
                    # ylims =(0.,Ylim)
                    )
        scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal")
        B2_plt = plot(sol, vars = [m2_HKCI, m2_PhlF],
                    label = ["Q2: output" L"$\bar{Q}$2: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The second bit counter",
                    # ylims =(0.,Ylim)
                    )
        B3_plt = plot(sol, vars = [m3_HKCI, m3_PhlF],
                    label = ["Q3: output" L"$\bar{Q}$3: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The third bit counter",
                    # ylims =(0.,Ylim)
                    )
        B4_plt = plot(sol, vars = [m4_HKCI, m4_PhlF],
                    label = ["Q4: output" L"$\bar{Q}$4: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The fourth bit counter",
                    # ylims =(0.,Ylim)
                    )
        B5_plt = plot(sol, vars = [m5_HKCI, m5_PhlF],
                    label = ["Q5: output" L"$\bar{Q}$5: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The fifth bit counter",
                    # ylims =(0.,Ylim)
                    )
        Bit5_plt_B = plot(
                        B1_plt, B2_plt, B3_plt, B4_plt, B5_plt,
                        layout = (5,1),
                        xtickfontsize=11, ytickfontsize=11, size=(1000,800),
                        legend = :topright,
                        # legend = :outertopright,
                         legendfont = font("Times new roman", 10),
                        xlabel = "Time", ylabel = "Concentration", dpi = 400,
                        )
        return Bit5_plt_B
    end
end

## run some ploting
function single_run(;δ = 270, n =2.81, show = "yes")
    sol, ts = run_prob_5bits(;init_relax = 20000.,duration= δ , relax=20000., signal=20.0, K=0.081, n= n, up = 1.5, cycle = 40)
    # plot(sol, vars = g23)
    if show == "yes"
        plt_B = plot_5bit(sol,ts,up, opt = 'B')
        plt_C = plot_5bit(sol,ts,up, opt = 'C')
        display(plot(plt_B,plt_C, layout = (1,2)))
    else show == "no"
        return sol, ts
    end
end

single_run(δ = 260, n= 2.7,show = "yes") # works for 5 bits


for δ ∈ 200:5:350
    @show δ
    single_run(δ = δ, n= 2.7,show = "yes")
end





B1_plt = plot(sol, vars = [m1_HKCI, m1_PhlF],
                    label = ["Q1: output" L"$\bar{Q}$1: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The first bit counter",
                    # ylims =(0.,Ylim)
                    )
scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal")
B2_plt = plot(sol, vars = [m2_HKCI, m2_PhlF],
            label = ["Q2: output" L"$\bar{Q}$2: variable to feed back"],
            legend = :topright,
            legendtitle = "The second bit counter",
            # ylims =(0.,Ylim)
            )
B3_plt = plot(sol, vars = [m3_HKCI, m3_PhlF],
            label = ["Q3: output" L"$\bar{Q}$3: variable to feed back"],
            legend = :topright,
            legendtitle = "The third bit counter",
            # ylims =(0.,Ylim)
            )
B4_plt = plot(sol, vars = [m4_HKCI, m4_PhlF],
            label = ["Q4: output" L"$\bar{Q}$4: variable to feed back"],
            legend = :topright,
            legendtitle = "The fourth bit counter",
            # ylims =(0.,Ylim)
            )

B5_plt = plot(sol, vars = [m5_HKCI, m5_PhlF],
                    label = ["Q5: output" L"$\bar{Q}$5: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The fifth bit counter",
                    # ylims =(0.,Ylim)
                    )
B6_plt = plot(sol, vars = [m6_HKCI, m6_PhlF],
                    label = ["Q6: output" L"$\bar{Q}$6: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The sixth bit counter",
                    # ylims =(0.,Ylim)
                    )
B7_plt = plot(sol, vars = [m7_HKCI, m7_PhlF],
                    label = ["Q7: output" L"$\bar{Q}$7: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The seventh bit counter",
                    # ylims =(0.,Ylim)
                    )
B8_plt = plot(sol, vars = [m8_HKCI, m8_PhlF],
                    label = ["Q8: output" L"$\bar{Q}$8: variable to feed back"],
                    legend = :topright,
                    legendtitle = "The eigth bit counter",
                    # ylims =(0.,Ylim)
                    )