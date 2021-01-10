#Author: Tianchi Chen

using DifferentialEquations, ParameterizedFunctions,Plots;plotly()
using Latexify
using CSV, DataFrames

mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end

t1 = 1000.; t2 = 1040;
const tstop1 = [t1]; const tstop2 = [t2]; const tstop3 = [2000.]; const tstop4 = [2500.]
# condition functions
function condition(u,t,integrator)
  t in tstop1
end
function condition2(u,t,integrator)
  t in tstop2
end
# function condition3(u,t,integrator)
#   t in tstop3
# end
# function condition4(u,t,integrator)
#   t in tstop4
# end

function affect!(integrator)
  for c in full_cache(integrator)
    c.f1 = 20.5
  end
end
function affect2!(integrator)
  for c in full_cache(integrator)
    c.f1 = 0.0
  end
end
# function affect3!(integrator)
#   for c in full_cache(integrator)
#     c.f1 = 0.5
#   end
# end
# function affect4!(integrator)
#   for c in full_cache(integrator)
#     c.f1 = 0.0
#   end
# end




save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions);
cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions);
# cb3 = DiscreteCallback(condition3, affect3!, save_positions=save_positions);
# cb4 = DiscreteCallback(condition4, affect4!, save_positions=save_positions);
# cbs = CallbackSet(cb,cb2,cb3,cb4)
cbs = CallbackSet(cb,cb2)


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

# Define functions
γ = 0.025
ξ = 0.025



function plasmid1(du,u,p,t)
    #  H1_HlyIIR   ---> NOR (pBAD, BetI)
    du[1] = ξ*(HlyIIR_min + (HlyIIR_max - HlyIIR_min)*HlyIIR_K^HlyIIR_n/(HlyIIR_K^HlyIIR_n + (u[7] + u.f1)^HlyIIR_n)) - γ*(u[1])

    #  N1_LmrA    ----> NOR (pBAD, HlyIIR)
    du[2] = ξ*(LmrA_min + (LmrA_max - LmrA_min)*LmrA_K^LmrA_n/(LmrA_K^LmrA_n + (u[1] + u.f1)^LmrA_n)) - γ*(u[2])

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

u0 = SimType([0.0;0.0;0.0;0.0;0.0;1.0;0.0], 0.0)
prob = ODEProblem(plasmid1,u0,(0.0,2900.0))
tstop = [t1;t2]# const
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
plot(sol, vars=[(0,6),(0,7)], xlabel = "time", ylabel = "concentration")#,vars=(6)

# control parameters
[sol[i].f1 for i in 1:length(sol)]
plot!(sol.t, [sol[i].f1 for i in 1:length(sol)], m = :cross)








































# Model ------------ without inputs

const t1 = [500]; const t2 = [550]; #const tstop3 = [15.]

mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::T
end

ts = [500,550, 1000,1200]
condition(u,t,integrator) = t in ts
function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p = 0.5
  elseif integrator.t == ts[2]
      integrator.p = 0.0
  elseif integrator.t == ts[3]
      integrator.p = 0.5
  elseif integrator.t == ts[4]
      integrator.p = 0.0
  end
end
save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions);


condition(u,t,integrator) = t in t1
condition2(u,t,integrator) = t in t2
# condition3(u,t,integrator) = t in tstop3
function affect!(integrator)
  integrator.p = 100.0
end
function affect2!(integrator)
  integrator.p = 0.0
end
# function affect3!(integrator)
#   integrator.p = 100.0
# end
save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions); cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions);# cb3 = DiscreteCallback(condition3, affect3!, save_positions=save_positions)
# cbs = CallbackSet(cb,cb2,cb3)
cbs = CallbackSet(cb,cb2)

γ = 0.025
ξ = 0.025
response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
degradation(x) = γ*x


SR = @ode_def_bare SR_latch begin
    # 1.
    dm_HlyIIR = ξ*response(HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n, m_BetI + p) - degradation(m_HlyIIR)
    # 2.
    dm_LmrA = ξ*response(LmrA_min, LmrA_max, LmrA_K, LmrA_n, m_HlyIIR + p) - degradation(m_LmrA)
    # 3.
    dm_SrpR = ξ*response(SrpR_min, SrpR_max, SrpR_K, SrpR_n, m_HlyIIR + m_BetI) - degradation(m_SrpR)
	# 4.
    dm_BM3R1 = ξ*response(BM3R1_min, BM3R1_max, BM3R1_K, BM3R1_n, m_LmrA + m_SrpR) - degradation(m_BM3R1)
    # 5.
    dm_PhIF = ξ*response(PhIF_min, PhIF_max, PhIF_K, PhIF_n, m_BM3R1) - degradation(m_PhIF)
    # 6.
    dm_PsrA = ξ*response(PsrA_min, PsrA_max, PsrA_K, PsrA_n, m_PhIF + m_BetI ) - degradation(m_PsrA)
    # 7.
    dm_BetI = ξ*response(BetI_min, BetI_max, BetI_K, BetI_n, m_BM3R1 + m_PsrA) - degradation(m_BetI)
end
u0 = SType([0.0; 0; 0; 0; 0; 0; 0], 0.0)
# u0 = [0.0; 0; 0; 0; 0; 0; 0]
p = 0.0
prob = ODEProblem(SR,u0,(0.0,1000.0),p)
sol = solve(prob,Tsit5(),callback = cb, tstops=ts)
# sol = solve(prob,Tsit5())
plot(sol,vars=[:m_HlyIIR,:m_BetI,:m_PsrA], xlabel = "time", ylabel = "concentration")
