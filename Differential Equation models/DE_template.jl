
# --------------------------------------      ParameterizedFunctions model     -------------------------------------
using DifferentialEquations, ParameterizedFunctions,Plots;plotly()
using Latexify
using CSV, DataFrames

const tstop1 = [5.]; const tstop2 = [10.]; const tstop3 = [15.]

mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::T
end

condition(u,t,integrator) = t in tstop1
condition2(u,t,integrator) = t in tstop2
condition3(u,t,integrator) = t in tstop3
function affect!(integrator)
  integrator.p = 100.0
end
function affect2!(integrator)
  integrator.p = 0.0
end
function affect3!(integrator)
  integrator.p = 100.0
end
save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions); cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions); cb3 = DiscreteCallback(condition3, affect3!, save_positions=save_positions)
cbs = CallbackSet(cb,cb2,cb3)



# Import gate Parameters
para_s1 = CSV.read("database/para_s1.csv");
para_s4 = CSV.read("database/para_s4.csv");
p1 = para_s1[:,[:repressor, :Y_min, :Y_max, :K, :n]];
p4 = para_s4[:,[:repressor, :Y_min, :Y_max, :K, :n]];

# Load parameters for each gate
HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n = convert(Matrix, p1[p1.repressor.=="HlyIIR",[:Y_min, :Y_max, :K, :n]]) 
PhIF_min, PhIF_max, PhIF_K, PhIF_n = convert(Matrix, p1[p1.repressor.=="PhIF",[:Y_min, :Y_max, :K, :n]])[3,:]
PsrA_min, PsrA_max, PsrA_K, PsrA_n = convert(Matrix, p4[p4.repressor.=="PsrA",[:Y_min, :Y_max, :K, :n]]) 
LmrA_min, LmrA_max, LmrA_K, LmrA_n = convert(Matrix, p4[p4.repressor.=="LmrA",[:Y_min, :Y_max, :K, :n]]) 
AmtR_min, AmtR_max, AmtR_K, AmtR_n = convert(Matrix, p1[p1.repressor.=="AmtR",[:Y_min, :Y_max, :K, :n]]) 


# Define functions
γ = 0.025
ξ = 0.025
response(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
degradation(x) = γ*x

# Model
SR = @ode_def_bare SR_latch begin
    # t =1
    #  $65: H1_HlyIIR ---> NOR (In: $62:PhIF, $64:PsrA)
    dm_HlyIIR = ξ*response(HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n, m_PhIF + m_PsrA ) - degradation(m_HlyIIR)
    #  $62: P3_PhIF   ----> NOR (In: $65:HlyIIR, $61:LmrA )
    dm_PhIF = ξ*response(PhIF_min, PhIF_max, PhIF_K, PhIF_n, m_HlyIIR + m_LmrA) - degradation(m_PhIF)
	#  $64: R1_PsrA   ----> NOR (In: $63:AmtR, $t)
    dm_PsrA = ξ*response(PsrA_min, PsrA_max, PsrA_K, PsrA_n, m_AmtR + p) - degradation(m_PsrA)
    #  $61: N1_LmrA ----> NOT (In: t)
    dm_LmrA = ξ*response(LmrA_min, LmrA_max, LmrA_K, LmrA_n, p) - degradation(m_LmrA)
    #  $63: A1_AmtR ----> NOT (In: $65:HlyIIR )
    dm_AmtR = ξ*response(AmtR_min, AmtR_max, AmtR_K, AmtR_n, m_HlyIIR) - degradation(m_AmtR)
end
u0 = SimType([0.0; 0; 0; 0; 0], 0.0)
p = 0.0
prob = ODEProblem(SR,u0,(0.0,30.0),p)
const tstop = [5.;10.;15.]
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
plot(sol,vars=[:m_HlyIIR], xlabel = "time", ylabel = "concentration")




# ------------------------------------  DEDataArrays Method ----------------------------------
# The DEDataArray{T} type allows one to add other "non-continuous" variables to an array, which can be useful in many modeling situations involving lots of events. To define an DEDataArray, make a type which subtypes DEDataArray{T} with a field x for the "array of continuous variables" for which you would like the differential equation to treat directly. The other fields are treated as "discrete variables"

using DifferentialEquations, ParameterizedFunctions,Plots;plotly()
using Latexify
using CSV, DataFrames

const tstop1 = [100.]; const tstop2 = [600.]; const tstop3 = [1000.]; const tstop4 = [1500.]
# condition functions
function condition(u,t,integrator)
  t in tstop1
end
function condition2(u,t,integrator)
  t in tstop2
end
function condition3(u,t,integrator)
  t in tstop3
end
function condition4(u,t,integrator)
  t in tstop4
end

function affect!(integrator)
  for c in full_cache(integrator)
    c.f1 = 5.0
  end
end
function affect2!(integrator)
  for c in full_cache(integrator)
    c.f1 = 0.0
  end
end
function affect3!(integrator)
  for c in full_cache(integrator)
    c.f1 = 5.0
  end
end
function affect4!(integrator)
  for c in full_cache(integrator)
    c.f1 = 0.0
  end
end

mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end


save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions);cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)
cb3 = DiscreteCallback(condition3, affect3!, save_positions=save_positions);cb4 = DiscreteCallback(condition4, affect4!, save_positions=save_positions)
cbs = CallbackSet(cb,cb2,cb3,cb4)

# Import gate Parameters
para_s1 = CSV.read("database/para_s1.csv");
para_s4 = CSV.read("database/para_s4.csv");
p1 = para_s1[:,[:repressor, :Y_min, :Y_max, :K, :n]];
p4 = para_s4[:,[:repressor, :Y_min, :Y_max, :K, :n]];

# Load parameters for each gate
HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n = convert(Matrix, p1[p1.repressor.=="HlyIIR",[:Y_min, :Y_max, :K, :n]]) 
PhIF_min, PhIF_max, PhIF_K, PhIF_n = convert(Matrix, p1[p1.repressor.=="PhIF",[:Y_min, :Y_max, :K, :n]])[3,:]
PsrA_min, PsrA_max, PsrA_K, PsrA_n = convert(Matrix, p4[p4.repressor.=="PsrA",[:Y_min, :Y_max, :K, :n]]) 
LmrA_min, LmrA_max, LmrA_K, LmrA_n = convert(Matrix, p4[p4.repressor.=="LmrA",[:Y_min, :Y_max, :K, :n]]) 
AmtR_min, AmtR_max, AmtR_K, AmtR_n = convert(Matrix, p1[p1.repressor.=="AmtR",[:Y_min, :Y_max, :K, :n]]) 

γ = 0.025
ξ = 0.025


function plasmid1(du,u,p,t)
    #  $65: H1_HlyIIR ---> NOR (In: $62:PhIF, $64:PsrA)
    du[1] = ξ*(HlyIIR_min + (HlyIIR_max - HlyIIR_min)*HlyIIR_K^HlyIIR_n/(HlyIIR_K^HlyIIR_n + (u[2]+u[3])^HlyIIR_n)) - γ*(u[1])

    #  $62: P3_PhIF   ----> NOR (In: $65:HlyIIR, $61:LmrA )
    du[2] = ξ*(PhIF_min + (PhIF_max - PhIF_min)*PhIF_K^PhIF_n/(PhIF_K^PhIF_n + (u[1] + u[4])^PhIF_n)) - γ*(u[2])

	#  $64: R1_PsrA   ----> NOR (In: $63:AmtR, $t)
    du[3] = ξ*(PsrA_min + (PsrA_max - PsrA_min)*PsrA_K^PsrA_n/(PsrA_K^PsrA_n + (u[5] + u.f1)^PsrA_n)) - γ*(u[3])

    #  $61: N1_LmrA ----> NOT (In: u.f1)
    du[4] = ξ*(LmrA_min + (LmrA_max - LmrA_min)*LmrA_K^LmrA_n/(LmrA_K^LmrA_n + (u.f1)^LmrA_n)) - γ*(u[4])

    #  $63: A1_AmtR ----> NOT (In: $65:HlyIIR )
    du[5] = ξ*(AmtR_min + (AmtR_max - AmtR_min)*AmtR_K^AmtR_n/(AmtR_K^AmtR_n + (u[1])^AmtR_n)) - γ*(u[5])
end

u0 = SimType([0.0;0.0;0.0;0.0;0.0], 0.0)
prob = ODEProblem(plasmid1,u0,(0.0,3025.0))
const tstop = [100.;600.;1000.;1500.]
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
plot(sol,vars=(1), xlabel = "time", ylabel = "concentration")
# control parameters
[sol[i].f1 for i in 1:length(sol)]








#  ------------------------------------------------------------------------------------------------
#  --------------------    test the toggle switch part only $62: PhIF;  $65: HlyIIR  🔺 works! ---------------------
using DifferentialEquations, ParameterizedFunctions,Plots;plotly()
using Latexify
using CSV, DataFrames

# Import gate Parameters
para_s1 = CSV.read("database/para_s1.csv");
para_s4 = CSV.read("database/para_s4.csv");
p1 = para_s1[:,[:repressor, :Y_min, :Y_max, :K, :n]];
p4 = para_s4[:,[:repressor, :Y_min, :Y_max, :K, :n]];

# Load parameters for each gate
HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n = convert(Matrix, p1[p1.repressor.=="HlyIIR",[:Y_min, :Y_max, :K, :n]]) 
PhIF_min, PhIF_max, PhIF_K, PhIF_n = convert(Matrix, p1[p1.repressor.=="PhIF",[:Y_min, :Y_max, :K, :n]])[3,:]



mutable struct SimType2{T} <: DEDataVector{T}
    x::Array{T,1}
    uA::T
    uB::T
end
const tstop1 = [100.]; const tstop2 = [400.]; const tstop3 = [700.]

condition(u,t,integrator) = t in tstop1
condition2(u,t,integrator) = t in tstop2
condition3(u,t,integrator) = t in tstop3
function affect!(integrator)
  for c in full_cache(integrator)
  c.uA = 1.0
  c.uB = 0.0
end
end
function affect2!(integrator)
  for c in full_cache(integrator)
  c.uA = 0.0
  c.uB = 1.0
end
end
function affect3!(integrator)

  for c in full_cache(integrator)
  c.uA = 1.0
  c.uB = 0.0
end
end
save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions); cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions); cb3 = DiscreteCallback(condition3, affect3!, save_positions=save_positions)
cbs = CallbackSet(cb,cb2,cb3)

γ = 0.025
ξ = 0.025



function f(du,u,p,t)
    du[1] = ξ*(HlyIIR_min + (HlyIIR_max - HlyIIR_min)*HlyIIR_K^HlyIIR_n/(HlyIIR_K^HlyIIR_n + (u[2]+u.uB)^HlyIIR_n)) - γ*(u[1])
    du[2] = ξ*(PhIF_min + (PhIF_max - PhIF_min)*PhIF_K^PhIF_n/(PhIF_K^PhIF_n + (u[1] + u.uA)^PhIF_n)) - γ*(u[2])
end

u0 = SimType2([0.0;0.0], 0.0, 0.0)
prob = ODEProblem(f,u0,(0.0,1000.0))
const tstop = [100.;400.;700.]
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)

plot(sol)

[sol[i].uA for i in 1:length(sol)]
[sol[i].uB for i in 1:length(sol)]
