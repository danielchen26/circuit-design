## Package
using DataFrames, CSV
using DifferentialEquations, ModelingToolkit
using Plots; gr(fontfamily = "Souce Code Pro for Powerline");
include("functions.jl")
## ==== Build multiple counter connectors : 2 Bits counter case 📗 =========
### Define a differential equation system
## ======= Build multiple counter connectors : 2 Bits counter case 📗 =========
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
### Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t) g1(t) g2(t) g3(t) m2_LexA1(t) m2_IcaR(t) m2_CI1(t) m2_PsrA(t) m2_BM3RI(t) m2_HKCI(t) m2_PhlF(t)
@derivatives D'~t
##
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
##
de2 = ODESystem(eqs2, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF, g1, g2, g3, m2_LexA1, m2_IcaR, m2_CI1, m2_PsrA,m2_BM3RI, m2_HKCI, m2_PhlF], [up,dn,K,n,γ,ξ,p])
ode_f2 = ODEFunction(de2)
## Solve DifferentialEquations
u0 = rand(1:22., length(ode_f2.syms))
Δ0, Δ, δ, cycle, A, tspan, time, signal, ts, cb, p = init_control(Δ0 = 20000., Δ = 20000., δ = 270, cycle = 20)
param = [1.5,0.002,0.081,2.81,0.025,0.025,p]
prob0 = ODEProblem(ode_f2, u0, tspan, param)
sol = solve(prob0, Tsit5(), callback = cb, tstops = ts, reltol = 1e-13, abstol = 1e-16)

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
sol, ts = run_prob_2bits(;init_relax = 5000., duration=270.,relax=5000.,signal=20.,K=0.081,n=2.81, up= 1.5, cycle=20)

## Plots for paper
function plot_2bit(sol, ts, up, param)
    Ylim = 3*up
    C_plt = plot(sol, vars = [:g3], legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    B1_plt = plot(sol, vars = [:m1_HKCI, :m1_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    scatter!(B1_plt, ts, 2*ones(length(ts)),label="Signal")
    B2_plt = plot(sol, vars = [:m2_HKCI, :m2_PhlF],legend = :topright, legendfontsize	= 7, ylims =(0.,Ylim))
    Bit2_plt = plot(C_plt,B1_plt,B2_plt,layout = (3,1),xtickfontsize=15,ytickfontsize=15, size=(1000,800),legend = :topright, legendfont = font("Times new roman", 8), title = ["C1 & Param: $param" "B1" "B2"])
    display(Bit2_plt)
end
plot_2bit(sol, ts, 1.5, param)




## check B3 switching
## Building cost function for param
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
ts_ID, C1 = Carrying_bit(sol, ts)

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
##
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
##
function osci(G6, up)
    thred = up/2
    G6b = G6 .> thred
    G6b1 = [G6b[i] for i = 1:2:length(G6b)]
    G6b2 = [G6b[i] for i = 2:2:length(G6b)]
    # @show G6b1, G6b2, G6b1 .+ G6b2
    return sum(G6b1 .+ G6b2), sum(G6b1 .* G6b2) # comp should be 4, and sum should be 0
end

Δ0 = 5000.

function cost_bit2(sol, ts, up)
    thred = up/2
    # set T0 for multibit counter
    ts_ID, C1 = Carrying_bit(sol, ts)
    T0 = MB_T0(ts_ID, C1)
    @show T0
    if length(T0) > 0
        G6,  G7 = Switch_cost(sol, Δ0, ts, T0; BNO = "B1")
        G16, G17 = Switch_cost(sol, Δ0, ts, T0; BNO ="B2")
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



## Check 3Bits_DB
db2 = CSV.read("2Bits_DB.csv")
db2
for i =1:size(db2)[1]
    # i = rand(1:size(db2)[1])
    δδ = db2[i,:].δ; AA =db2[i,:].A;  KK = db2[i,:].K; nn = db2[i,:].n; upp = db2[i,:].up
    sol, ts = run_prob_2bits(;init_relax = 5000., duration=δδ,relax=5000.,signal=AA, K=KK, n=nn, up= upp, cycle=20)
    param = [KK, nn, δδ, AA, upp]
    plot_2bit(sol, ts, upp, param)
    cost_bit2(sol, ts, upp)
end
