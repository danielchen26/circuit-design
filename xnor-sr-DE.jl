using DifferentialEquations, ParameterizedFunctions,Plots;plotly()
using Latexify
using CSV, DataFrames

# # General parameters
# γ = 0.025
# ξ = 0.025
#
# #functions
# respone(min, max, K, x, n) = (min + (max - min)*K^n/(K^n + x^n))
# degradation(x) = γ*x
#
# # SR latch parameter
# min_x = 1
# max_x = 1
# Kx = 1
# nx = 1
#
# min_y = 1
# max_y = 1
# Ky = 1
# ny = 1
#
# # Model
# SR = @ode_def_bare SR_latch begin
#     dm_qa = respone(min_x, max_x, Kx, m_qb, nx)
#     dm_qb = respone(min_y, max_y, Ky, m_qa, ny)
#
# end
#
# u0 = [1.0, 0.0]
# tspan = (0.0,1.0)
# prob = ODEProblem(SR,u0,tspan)
# sol = solve(prob)
# plot(sol)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Load parameter Pandas version
# para = read_csv("para.csv");
# para0= read_csv("para0.csv");
# para0
# iloc(para0)[1,:]
# p0[p0[:repressor]=="AmeR"][[:Y_min, :Y_max,:K,:n]]


# Import gate Parameters
para_s1 = CSV.read("para_s1.csv");
para_s4 = CSV.read("para_s4.csv");
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
respone(min, max, K, n, x) = (min + (max - min)*K^n/(K^n + x^n))
degradation(x) = γ*x


# DifferentialEquations models for plasmid1

pl1 = @ode_def_bare plasmid1 begin
    # t =1
    #  $65: H1_HlyIIR ---> NOR (In: $62:PhIF, $64:PsrA)
    dm_HlyIIR = ξ*respone(HlyIIR_min, HlyIIR_max, HlyIIR_K, HlyIIR_n, m_PhIF + m_PsrA ) - degradation(m_HlyIIR)

    #  $62: P3_PhIF   ----> NOR (In: $65:HlyIIR, $61:LmrA )
    dm_PhIF = ξ*respone(PhIF_min, PhIF_max, PhIF_K, PhIF_n, m_HlyIIR + m_LmrA) - degradation(m_PhIF)

	#  $64: R1_PsrA   ----> NOR (In: $63:AmtR, $t)
    dm_PsrA = ξ*respone(PsrA_min, PsrA_max, PsrA_K, PsrA_n, m_AmtR + t) - degradation(m_PsrA)

    #  $61: N1_LmrA ----> NOT (In: t)
    dm_LmrA = ξ*respone(LmrA_min, LmrA_max, LmrA_K, LmrA_n, t) - degradation(m_LmrA)

    #  $63: A1_AmtR ----> NOT (In: $65:HlyIIR )
    dm_AmtR = ξ*respone(AmtR_min, AmtR_max, AmtR_K, AmtR_n, m_HlyIIR) - degradation(m_AmtR)
end



u0 = [1.0, 0.0, 0, 0, 0]
tspan = (0.0,1.0)
prob = ODEProblem(pl1,u0,tspan)
sol = solve(prob)
plot(sol)




latexify(pl1)
