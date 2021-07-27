## ============ The 7 circuits for 1bit counter with physical gates
   
#  ============= The replacement table ==============
# if to replace 1 cello gate
# │ 1   │ AmeR      │ F1     │ 0.2     │ 3.8     │ 0.09    │ 1.4     │
# │ 2   │ AmtR      │ A1     │ 0.06    │ 3.8     │ 0.07    │ 1.6     │
# │ 8   │ IcaRA     │ I1     │ 0.08    │ 2.2     │ 0.1     │ 1.4     │
# │ 9   │ LitR      │ l1     │ 0.07    │ 4.3     │ 0.05    │ 1.7     │
# │ 15  │ QacR      │ Q1     │ 0.01    │ 2.4     │ 0.05    │ 2.7     │
# │ 16  │ QacR      │ Q2     │ 0.03    │ 2.8     │ 0.21    │ 2.4     │

# 1. cello  gates in order: => [HlyIIR, LmrA, SrpR, BM3R1, PhIF, PsrA,BetI]
p_unique1 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
circuit1_name = ["H1HlyIIR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
# 2. cello_replace_[0.2, 3.8, 0.09, 1.4]_2   🍏 LmrA ==> AmeR
p_unique2 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 3.8, 0.09, 1.4],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
circuit2_name = ["H1HlyIIR", "F1AmeR", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
# 3. cello_replace_[0.01, 2.4, 0.05, 2.7]_5  🍏 PhIF ==> QacR
p_unique3 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 2.4, 0.05, 2.7],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
circuit3_name = ["H1HlyIIR", "N1LmrA", "S4SrpR", "B3BM3R1", "Q1QacR", "R1PsrA","E1BetI"]
# 4. 🔴 cello_replace_[0.03, 2.8, 0.21, 2.4]_1  🍏 HlyIIR ==> QacR
p_unique4 = Any[[0.03, 2.8, 0.21, 2.4],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
circuit4_name = ["Q2QacR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
# 5. cello_replace_[0.03, 2.8, 0.21, 2.4]_7  🍏 BetI ==> QacR
p_unique5 = Any[[0.07, 2.5, 0.19, 2.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.03, 2.8, 0.21, 2.4]]
circuit5_name = ["H1HlyIIR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","Q2QacR"]
# 6. cello_replace_[0.06, 3.8, 0.07, 1.6]_1  🍏 HlyIIR ==> AmtR
p_unique6 = Any[[0.06, 3.8, 0.07, 1.6],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
circuit6_name = ["A1AmtR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
# 7. cello_replace_[0.07, 4.3, 0.05, 1.7]_1  🍏 HlyIIR ==> LitR
p_unique7 = Any[[0.07, 4.3, 0.05, 1.7],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
circuit7_name = ["L1LitR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]



## 1. Use GRUtils to plot, but the rgb color does not match with HEX right now.
using GRUtils, Colors
fig1 = Figure()  
draw(fig1)
x = 1e-3:0.01:1e3
Hill(dn, up, K, n, x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
y = Hill.(p_unique[2]...,x)
# Making a line plot is as simple as this:
#define color 
cl = [61, 104, 61]./255
plot(x,y, linewidth = 10, linecolor=GRUtils.color(cl...))
draw(fig1)
xlog(true)
ylog(true)
xlabel("BetI")
# xticklabels([""])
# yticklabels([""])
background(nothing)
grid(false)
ylim(1e-2,1e3)






## 2. Use the GR() backend to plot, but current only have box boader, but not as modular and professional as the GRUtils.
using Colors, Plots;gr() # the following plot need gr() backend.
using LaTeXStrings
# 4th one of the 7 circuits.
p_unique = Any[[0.03, 2.8, 0.21, 2.4],[0.2, 2.2, 0.18, 2.1],[0.007, 2.1, 0.1, 2.8],[0.01, 0.8, 0.26, 3.4],[0.01, 3.9, 0.03, 4.0],[0.2, 5.9, 0.19, 1.8],[0.07, 3.8, 0.41, 2.4]]
x = 1e-3:0.01:100
Hill(dn, up, K, n, x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
best_circuit_gates_name = ["Q2QacR", "N1LmrA", "S4SrpR", "B3BM3R1", "P1PhIF", "R1PsrA","E1BetI"]
best_circuit_gates_color =[colorant"#a0c54d", colorant"#d04d2b", colorant"#3d683e", colorant"#cd3330", colorant"#e5a237", colorant"#a94b97", colorant"#55559f"]
function plot_activation_func(gate, gate_name, color)
    plt = plot(x, Hill.(gate...,x),lw=20, 
			   color = color, 
			   xaxis=:log, yaxis=:log,
			   thickness_scaling = 1.3,
			#    bordercolor="white",
			   framestyle = :box,
			   grid=false,
			   minorgrid=true,
               minorticks = 10,
			#    labels = L"Hill(x) = dn + \frac{(up - dn) * K^{n}}{(K^{n} + x^{n})}",
               ylims =[1e-3,1e2],tickfont = Plots.font("Helvetica", 6), 
               guidefont = Plots.font("Helvetica", 16), 
               legendfontsize = 30, 
               background_color_legend = nothing,
            #    legend = false, 
               fg_legend = :transparent,
               label = gate_name,
               dpi = 300
              )
    display(plt)
    return plt
end 
plot_activation_func(p_unique[1], best_circuit_gates_name[1], best_circuit_gates_color[1])
# to save these 7 activation fucntions to Paper_plots/physical_gate_cases/
for i in 1:7
    gate = p_unique[i]
    name = best_circuit_gates_name[i]
    color = best_circuit_gates_color[i]
    plt = plot_activation_func(gate, name, color)
    savepath = "./DEmodels/scripts/Paper_plots/physical_gate_cases/best_circuit_hill/Gate_$i"*"_name_:"*"$name"*".svg"
    savefig(plt,savepath)
end


# Save all 7 circuits' activation functions listed above in 7 folders
using Colors, Plots;gr() # the following plot need gr() backend.
using LaTeXStrings
circuit_set = [p_unique1, p_unique2, p_unique3, p_unique4,p_unique5, p_unique6, p_unique7]
circuit_name = [circuit1_name, circuit2_name, circuit3_name, circuit4_name, circuit5_name, circuit6_name, circuit7_name]
circuit_gates_color =[colorant"#a0c54d", colorant"#d04d2b", colorant"#3d683e", colorant"#cd3330", colorant"#e5a237", colorant"#a94b97", colorant"#55559f"]
x = 1e-3:0.01:100
Hill(dn, up, K, n, x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
for j in 1:7
    circuit = circuit_set[j]
    j_set_name = circuit_name[j]
    dir_name = joinpath(pwd()*"/DEmodels/scripts/Paper_plots/physical_gate_cases/circuit_$j")
    if !ispath(dir_name)
        mkpath(dir_name)
        @show dir_name
    end
    for i in 1:7
        gate = circuit[i]
        name = j_set_name[i]
        color = circuit_gates_color[i]
        plt = plot_activation_func(gate, name, color)
        # savepath = "./DEmodels/scripts/Paper_plots/physical_gate_cases/best_circuit_hill/Gate_$i"*"_name_:"*"$name"*".svg"
        savepath = joinpath(dir_name*"/Gate_$i"*"_name_:"*"$name"*".png")
        savefig(plt,savepath)
    end
end 