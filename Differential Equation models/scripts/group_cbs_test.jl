
using DifferentialEquations, ParameterizedFunctions,Plots;pgfplots()
using Latexify
using CSV, DataFrames


mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::T
end



ts = [1500,2000]
condition(u,t,integrator) = t in ts
function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p = 0.0
  elseif integrator.t == ts[2]
      integrator.p = 0.0
  # elseif integrator.t == ts[3]
  #     integrator.p = 1.5
  # elseif integrator.t == ts[4]
  #     integrator.p = 0.0
  end
end
save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions);

# Import gate Parameters
para = CSV.read("database/para_s4.csv");
p4 = para[:,[:repressor, :Y_min, :Y_max, :K, :n]];



# Load parameters for each gate
HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n = convert(Matrix, p4[p4.repressor.=="HlyIIR",[:Y_min, :Y_max, :K, :n]]) 
PhIF_min, PhIF_max, PhIF_K, PhIF_n = convert(Matrix, p4[p4.repressor.=="PhIF",[:Y_min, :Y_max, :K, :n]])[1,:]
PsrA_min, PsrA_max, PsrA_K, PsrA_n = convert(Matrix, p4[p4.repressor.=="PsrA",[:Y_min, :Y_max, :K, :n]]) 
LmrA_min, LmrA_max, LmrA_K, LmrA_n = convert(Matrix, p4[p4.repressor.=="LmrA",[:Y_min, :Y_max, :K, :n]]) 
SrpR_min, SrpR_max, SrpR_K, SrpR_n = convert(Matrix, p4[p4.repressor.=="SrpR",[:Y_min, :Y_max, :K, :n]])[4,:]
BM3R1_min, BM3R1_max, BM3R1_K, BM3R1_n = convert(Matrix, p4[p4.repressor.=="BM3R1",[:Y_min, :Y_max, :K, :n]])[3,:]
BetI_min, BetI_max, BetI_K, BetI_n = convert(Matrix, p4[p4.repressor.=="BetI",[:Y_min, :Y_max, :K, :n]])

γ = 0.025
ξ = 0.025
# response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
# degradation(x) = γ*x

function plasmid1(du,u,p,t)
    #  H1_HlyIIR   ---> NOR (pBAD, BetI)
    du[1] = ξ*(HlyIIR_min + (HlyIIR_max - HlyIIR_min)*HlyIIR_K^HlyIIR_n/(HlyIIR_K^HlyIIR_n + (u[7] + p)^HlyIIR_n)) - γ*(u[1])
    #  N1_LmrA    ----> NOR (pBAD, HlyIIR)
    du[2] = ξ*(LmrA_min + (LmrA_max - LmrA_min)*LmrA_K^LmrA_n/(LmrA_K^LmrA_n + (u[1] + p)^LmrA_n)) - γ*(u[2])
	#  R1_SrpR    ----> NOR (HlyIIR, BetI)
    du[3] = ξ*(SrpR_min + (SrpR_max - SrpR_min)*SrpR_K^SrpR_n/(SrpR_K^SrpR_n + (u[1] + u[7])^SrpR_n)) - γ*(u[3])
    #  B3_BM3R1   ----> NOR (LmrA, SrpR)
    du[4] = ξ*(BM3R1_min + (BM3R1_max - BM3R1_min)*BM3R1_K^BM3R1_n/(BM3R1_K^BM3R1_n + (u[2] + u[3])^BM3R1_n)) - γ*(u[4])
    #  P3_PhIF    ----> NOT (BM3R1 )
    du[5] = ξ*(PhIF_min + (PhIF_max - PhIF_min)*PhIF_K^PhIF_n/(PhIF_K^PhIF_n + (u[4])^PhIF_n)) - γ*(u[5])
    #  R1_PsrA    ----> NOR (PhIF, BetI)
    du[6] = ξ*(PsrA_min + (PsrA_max - PsrA_min)*PsrA_K^PsrA_n/(PsrA_K^PsrA_n + (u[5] + u[7])^PsrA_n)) - γ*(u[6])
    #  E1_BetI    ----> NOR (PsrA, BM3R1 )
    du[7] = ξ*(BetI_min + (BetI_max - BetI_min)*BetI_K^BetI_n/(BetI_K^BetI_n + (u[4] + u[6])^BetI_n)) - γ*(u[7])
end
# u0 = SType([3.0;0.0;1.0;0.4;1.0;0.0;1.0], 0.0)

u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
p=0.0
prob = ODEProblem(plasmid1,u0,(0.0,2200.0),p)
sol = solve(prob,Tsit5())#callback = cb, tstops=ts
plot(sol,vars=[(0,6),(0,7)], xlabel = "time", ylabel = "concentration")


u0_1 = SType(sol[end].x, 20.0)
p=20.0
prob2 = ODEProblem(plasmid1,u0_1,(0.0,2200.0),p)
sol2 = solve(prob2,Tsit5())#callback = cb, tstops=ts
plot(sol2,vars=[(0,6),(0,7)], lw =1,xlabel = "time", ylabel = "concentration")



for i in collect(1:10)
    u0 = SType(Float64[i for i in rand(1:22,7)], 0.0)
    p=0.0
    prob = ODEProblem(plasmid1,u0,(0.0,2200.0),p)
    sol = solve(prob,Tsit5(),callback = cb, tstops=ts)
    # display(plot!(sol,vars=[(0,6),(0,7)], xlabel = "time", ylabel = "concentration"))
    pp =  plot!(p, xlabel = "time", ylabel = "concentration",yticks = 0:0.5:10,)
end


# -----------note can not use ParameterizedFunctions
# SR = @ode_def_bare SR_latch begin
#     # 1.
#     dm_HlyIIR = ξ*response(HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n, m_BetI + p) - degradation(m_HlyIIR)
#     # 2.
#     dm_LmrA = ξ*response(LmrA_min, LmrA_max, LmrA_K, LmrA_n, m_HlyIIR + p) - degradation(m_LmrA)
#     # 3.
#     dm_SrpR = ξ*response(SrpR_min, SrpR_max, SrpR_K, SrpR_n, m_HlyIIR + m_BetI) - degradation(m_SrpR)
# 	# 4.
#     dm_BM3R1 = ξ*response(BM3R1_min, BM3R1_max, BM3R1_K, BM3R1_n, m_LmrA + m_SrpR) - degradation(m_BM3R1)
#     # 5.
#     dm_PhIF = ξ*response(PhIF_min, PhIF_max, PhIF_K, PhIF_n, m_BM3R1) - degradation(m_PhIF)
#     # 6.
#     dm_PsrA = ξ*response(PsrA_min, PsrA_max, PsrA_K, PsrA_n, m_PhIF + m_BetI ) - degradation(m_PsrA)
#     # 7.
#     dm_BetI = ξ*response(BetI_min, BetI_max, BetI_K, BetI_n, m_BM3R1 + m_PsrA) - degradation(m_BetI)
# end
# u0 = SType([0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0], 0.0)
# p = 0.0
# prob = ODEProblem(SR, u0, (0.0, 2000.0), p)
# sol = solve(prob,Tsit5(),callback = cb, tstops=ts)
# plot(sol, xlabel = "time", ylabel = "concentration")#,vars=[:m_HlyIIR,:m_BetI,:m_PsrA]
# # ploting response
# function t2p(t, ts)
#     idx = findfirst(x->x>t, @view(ts[:, 1])) # find the first idx where `ts[idx] > t`
#     idx === nothing && return ts[end, 2] # if there is no such point, we should use the last point
#     idx -= 1 # look at the previous point
#     iszero(idx) && return 0. # if `idx` is zero, we should use the initial condition
#     return ts[idx, 2]
# end
# x = collect(0:2:2000)
# plot!(x, t2p.(x, Ref(ts)))


# # working example (simple test)-----------------------------------------------
# using DifferentialEquations, ParameterizedFunctions,Plots;plotly()
# mutable struct SType{T} <: DEDataVector{T}
#     x::Array{T,1}
#     p::T
# end
#
# ts = [5,6, 10,11]
# condition(u,t,integrator) = t in ts
# function affect!(integrator)
#   if integrator.t == ts[1]
#       integrator.p = 3.0
#   elseif integrator.t == ts[2]
#       integrator.p = 0.0
#   elseif integrator.t == ts[3]
#       integrator.p = 3.0
#   elseif integrator.t == ts[4]
#       integrator.p = 0.0
#   end
# end
#
# save_positions = (true,true)
# cb = DiscreteCallback(condition, affect!, save_positions=save_positions);
#
# # -- ⃝ 1
# function f(du,u,p,t)
#     du[1] = -0.5*u[1] + p
#     du[2] = -0.5*u[2]
# end
# u0 = SType([10.0;10.0], 0.0)
# p = 0.0
# prob = ODEProblem(f,u0,(0.0,15.0),p)
# sol = solve(prob,Tsit5(),callback = cb, tstops=ts)
# plot(sol)
#
#
#
# # -- ⃝ 2
# test = @ode_def_bare test_model begin
#     dx = -0.5x + p
#     dv = -0.5v
# end p
#
# u0 = SType([10.0;10.0], 0.0)
# tspan = (0.0,15.0)
# p = 0
# prob = ODEProblem(test,u0,tspan,p)
# sol = solve(prob,Tsit5(),callback=cb,tstops=ts)
# plot(sol)
# # working example-----------------------------------------------



# # yingbo example --------
# using DifferentialEquations, ParameterizedFunctions, Plots; plotly()
# mutable struct SType{T} <: DEDataVector{T}
#     x::Array{T,1}
#     p::T
# end
#
# ts = [5 0.5
#       10 0.
#       20 0.5
#       25 0.]
# condition(u,t,integrator) = t in @view(ts[:, 1])
# function affect!(integrator)
#     idx = searchsortedfirst(@view(ts[:, 1]), integrator.t)
#     integrator.p = ts[idx, 2]
# end
# save_positions = (true,true)
# cb = DiscreteCallback(condition, affect!, save_positions=save_positions);
#
# function f(du,u,p,t)
#     du[1] = -0.5*u[1] + p
#     du[2] = -0.5*u[2]
# end
# u0 = SType([10.0;10.0], 0.0)
# p = 0.0
# prob = ODEProblem(f,u0,(0.0,30.0),p)
# sol = solve(prob,Tsit5(),callback = cb, tstops=ts)
# plot(sol)
#
#
#
#
# # yingbo functions to plot inpulse
# ts = [500 100
#       550 0
#       1000 100
#       1200 100]
# function t2p(t, ts)
#     idx = findfirst(x->x>t, @view(ts[:, 1])) # find the first idx where `ts[idx] > t`
#     idx === nothing && return ts[end, 2] # if there is no such point, we should use the last point
#     idx -= 1 # look at the previous point
#     iszero(idx) && return 0. # if `idx` is zero, we should use the initial condition
#     return ts[idx, 2]
# end
# x = collect(0:2000)
# plot!(x, t2p.(x, Ref(ts)), m = :dot)



# # ------yingbo's way can not replicat the graph from doc
# using DifferentialEquations, ParameterizedFunctions, Plots; plotly()
# mutable struct SType{T} <: DEDataVector{T}
#     x::Array{T,1}
#     p::T
# end
#
# ts = [5 1.5
#       8 -1.5]
# condition(u,t,integrator) = t in @view(ts[:, 1])
# function affect!(integrator)
#     idx = searchsortedfirst(@view(ts[:, 1]), integrator.t)
#     integrator.u.p = ts[idx, 2]
# end
# save_positions = (true,true)
# cb = DiscreteCallback(condition, affect!, save_positions=save_positions);
#
# function f(du,u,p,t)
#     du[1] = -0.5*u[1] + u.p
#     du[2] = -0.5*u[2]
# end
# u0 = SType([10.0;10.0], 0.0)
# prob = ODEProblem(f,u0,(0.0,10.0))
# sol = solve(prob,Tsit5(),callback = cb, tstops=ts)
# plot(sol)
