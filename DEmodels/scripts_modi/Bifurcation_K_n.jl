#Author: Tianchi Chen
## ------ Import package and functions
using DifferentialEquations, ModelingToolkit
using Plots
using Latexify, Random, Base
using CSV, DataFrames, ProgressMeter
using LaTeXStrings
using Peaks
using StatsPlots
# gr(fontfamily = "Souce Code Pro for Powerline");
pyplot()
# include("functions.jl")# using BlackBoxOptim, LinearAlgebra



## ==== Build multiple counter connectors : 1Bit counter case 📗 =========
hill(x) = dn + (up - dn) * K^n / (K^n + abs(x)^n)
deg(x) = γ * x
 # Define a differential equation system
@parameters t up dn K n γ ξ p
@variables m1_LexA1(t) m1_IcaR(t) m1_CI1(t) m1_PsrA(t) m1_BM3RI(t) m1_HKCI(t) m1_PhlF(t)
D = Differential(t)
eqs1 = [
    # Bit 1 =================
    D(m1_LexA1) ~ ξ * hill(m1_PhlF + p)        - deg(m1_LexA1),
    D(m1_IcaR ) ~ ξ * hill(m1_LexA1 + p)       - deg(m1_IcaR),
    D(m1_CI1  ) ~ ξ * hill(m1_LexA1 + m1_PhlF) - deg(m1_CI1),
    D(m1_PsrA ) ~ ξ * hill(m1_IcaR + m1_CI1)   - deg(m1_PsrA),
    D(m1_BM3RI) ~ ξ * hill(m1_PsrA)            - deg(m1_BM3RI),
    D(m1_HKCI ) ~ ξ * hill(m1_BM3RI + m1_PhlF) - deg(m1_HKCI),
    D(m1_PhlF ) ~ ξ * hill(m1_PsrA + m1_HKCI)  - deg(m1_PhlF)]
@named de1 = ODESystem(eqs1, t, [m1_LexA1, m1_IcaR, m1_CI1, m1_PsrA, m1_BM3RI, m1_HKCI,m1_PhlF], [up,dn,K,n,γ,ξ,p])
# ode_f1 = ODEFunction(de1)
##

function check_periodic(sol, idx)
    # pks, vals = findmaxima(sol[idx,:])
    pks_max, vals_max = findmaxima(sol[idx,:])
    pks_min, vals_min = findminima(sol[idx,:])
    # len = min(length(vals_min), length(vals_max))
    # vals_min = vals_min[1:len]; vals_max = vals_max[1:len]

    ref_value = vals_max[end] 
    ref_id = pks_max[end]
    high_diff = vals_max .- ref_value
    low_diff = vals_min .- ref_value
    osci = round(abs.(low_diff)[end], digits = 3)
    return pks_min, vals_min, pks_max, vals_max, osci
end









u0 = rand(7)
tspan = (0.0,1e4)
param = [1.5, 0.002, 0.051, 2.0, 0.025, 0.025, 20]
prob = ODEProblem(de1, u0, tspan, param)
sol0 = solve(prob)
plot(sol0, vars = [m1_HKCI,m1_PhlF])
pks_min, vals_min, pks_max, vals_max = check_periodic(sol0, 7)
scatter!(sol0.t[pks_min], vals_min)
scatter!(sol0.t[pks_max], vals_max)
##




##
using ProgressMeter
df_bifur = DataFrame(K = Float64[], n = Float64[], cri = [])
df_bifur_3d = DataFrame(K = Float64[], n = Float64[], up = Float64[], cri = [])
@showprogress for K = 0.011: 0.02: 1.0, n = 1.5:0.1:10.0, up = 1:10
    u0 = rand(7)
    tspan = (0.0,1e4)
    param = [up, 0.002, K, n, 0.025, 0.025, 20] # for 2d up set to 1.5
    prob = ODEProblem(de1, u0, tspan, param)
    sol0 = solve(prob)

    # check Q 
    # pks_min_6, vals_min_6, pks_max_6, vals_max_6 = check_periodic(sol0, 6)
    # # make the high oscillation values as a reference point
    # ref_value_6 = vals_max_6[end] 
    # ref_id_6 = pks_max_6[end]
    # high_diff_6 = vals_max_6 .- ref_value_6
    # low_diff_6 = vals_min_6 .- ref_value_6
    # osci_6 = round(abs.(low_diff_6)[end], digits = 3)
    # @show osci_6
    _, _, _, _, osci_6= check_periodic(sol0, 6)
    @show osci_6
    
    # check \bar{Q }
    # pks_min_7, vals_min_7, pks_max_7, vals_max_7 = check_periodic(sol0, 7)
    # # make the high oscillation values as a reference point
    # ref_value_7 = vals_max_7[end] 
    # ref_id_7 = pks_max_7[end]
    # high_diff_7 = vals_max_7 .- ref_value_7
    # low_diff_7 = vals_min_7 .- ref_value_7
    # osci_7 = round(abs.(low_diff_7)[end], digits = 3)
    # @show osci_7

    _, _, _, _, osci_7= check_periodic(sol0, 7)
    @show osci_7

    cri = osci_6 > 0.01 || osci_7 > 0.01
    @show cri 

    # for 2d plots
    # record = [K, n, cri]
    # cri == true ? push!(df_bifur, record) : push!(df_bifur, [K,n, 0.0])
    # for 3d plots
    record = [K, n, up, cri]
    cri == true ? push!(df_bifur_3d, record) : push!(df_bifur_3d, [K,n, up, 0.0])
    


    # plots 
    # plt = plot(sol0,vars = [m1_HKCI,m1_PhlF],
    #            legend = :outertopright,
    #             title = "Oscillation: $cri, G6: $osci_6, G7: $osci_7")
    # scatter!( plt, [sol0.t[ref_id_6]], [ref_value_6], c = :orange, marker = (:star, 10))
    # scatter!( plt, [sol0.t[ref_id_7]], [ref_value_7], c = :purple, marker = (:square, 10))
    # plt1 = plot(high_diff_6, c= :orange, label = "Hight Q residue")
    # plt2 = plot(low_diff_6,  c = :orange, label = "Low Qb residue")
    # plt3 = plot(high_diff_7, c = :purple, label = "High Q residue")
    # plt4 = plot(low_diff_7,  c = :purple, label = "Low Qb residue")
    # display(plot(plt, plt1, plt2, plt3, plt4, layout = @layout([a; b c; d e])))
end 
## 


## ===== if not-oscillating, check with 10 iterations of random initial conditions =====
function check_non_oscillating(df)
    df_sub = filter(row -> row.cri == 0, df)
    @show df_sub
    df_update = DataFrame(K = Float64[], n = Float64[], up = Float64[], cri = [])
    for row in eachrow(df_sub)
        for i = 1:10
                u0 = rand(7)
                tspan = (0.0,1e5)
                param = [row.up, 0.002, row.K, row.n, 0.025, 0.025, 20] # for 2d up set to 1.5
                prob = ODEProblem(de1, u0, tspan, param)
                sol0 = solve(prob)
                _, _, _, _, osci_6= check_periodic(sol0, 6)
                # @show osci_6
                _, _, _, _, osci_7= check_periodic(sol0, 7)
                # @show osci_7
                cri = osci_6 > 0.01 || osci_7 > 0.01
                @show cri
                if cri == true
                    record = [row.K, row.n, row.up, cri]
                    push!(df_bifur_3d, record) 
                    plt = plot(sol0,vars = [m1_HKCI,m1_PhlF],legend = :outertopright,
                                 title = "Oscillation: $cri, G6: $osci_6, G7: $osci_7")
                    display(plt)
                    # scatter!( plt, [sol0.t[ref_id_6]], [ref_value_6], c = :orange, marker = (:star, 10))
                    # scatter!( plt, [sol0.t[ref_id_7]], [ref_value_7], c = :purple, marker = (:square, 10))
                    # plt1 = plot(high_diff_6, c= :orange, label = "Hight Q residue")
                    # plt2 = plot(low_diff_6,  c = :orange, label = "Low Qb residue")
                    # plt3 = plot(high_diff_7, c = :purple, label = "High Q residue")
                    # plt4 = plot(low_diff_7,  c = :purple, label = "Low Qb residue")
                    # display(plot(plt, plt1, plt2, plt3, plt4, layout = @layout([a; b c; d e])))
                    break
                else cri == false
                    continue
                end
        end
    end 
    return df_update
end 
df_update = check_non_oscillating(df_bifur_3d)




df_bifur
CSV.write("./DEmodels/scripts_modi/3d_bifur.csv",df_bifur_3d)



plt_bifur = @df df_bifur scatter(
                    :n,
                    :K,
                    xlabel = "n",
                    ylabel = "K",
                    label = ["Oscillation" "Non Oscillation"],
                    markersize = 4.2,
                    group =:cri,
                    dpi = 400)


plt_bifur_3d = @df df_bifur_3d scatter(
                    :n,
                    :K,
                    :up,
                    xlabel = "n",
                    ylabel = "K",
                    zlabel = "up",
                    label = ["Non Oscillation" "Oscillation"],
                    markersize = 2.2,
                    group =:cri,
                    dpi = 400)
# savefig(plt_bifur, "./DEmodels/scripts_modi/Vis/bifur.png")

using VegaLite
df_bifur |>
    @vlplot(
        :circle,
        x = :n,
        y = :K,
        color= :cri
    )

# =============== using Makie to plot 3d scatter plot ============================================================
# convert df to Tensor
using GLMakie
# load 3d bifur data
df_bifur_3d  = CSV.File("./DEmodels/scripts_modi/3d_bifur.csv") |> DataFrame
function df_to_tensor(df, Kr, nr, up_r)
    df_i = filter([:K, :n, :up] => (K, n, up) -> K == Kr && n == nr && up == up_r, df)
    # @show df_i
    return df_i.cri[1]
end 

K_r = 0.011: 0.02: 1.0; n_r = 1.5:0.1:10.0;up_r = 1:10

cube = [df_to_tensor(df_bifur_3d, K, n, up) for K = K_r, n =n_r , up = up_r]
cmap = RGBAf0.(to_colormap(:viridis), 1.0)
scene = volume(cube, algorithm = :absorption,transparency=true,colormap = cmap,
figure=(backgroundcolor=:white,) )


## make example
using GLMakie, ForwardDiff
let
    f(x,y) = -5*x*y*exp(-x^2-y^2)
    x = -1:0.05:1.0
    y = -1:0.05:1.0
    z = [f(i,j) for i in x, j in y]
    # This is the same function as above, just modified so that it will
    # work with ForwardDiff
    g(x,y) = [-5*x*y*exp(-x^2-y^2)]
    J(xx,yy) = ForwardDiff.jacobian(x -> g(x[1], x[2]), [xx,yy])
    field(i,j) = Point(J(i,j)[1], J(i,j)[2])

    zmin, zmax = minimum(z), maximum(z)
    cmap = :viridis
    function plot()
        fig = Figure(resolution = (1200,600))
        ax1 = Axis3(fig, aspect = (1,1,1), perspectiveness = 0.5,
            elevation = π/3.5, azimuth = 0.1π, )
        ax2 = Axis(fig, aspect = DataAspect(), xlabel = "x", ylabel = "y")
        surface!(ax1, x, y, z, colormap = cmap, colorrange = (zmin, zmax))
        contour3d!(ax1, x, y, z .+ 0.005, levels = 15, linewidth = 2, color=(:white, 0.5))
        wireframe!(ax1, x, y, z, overdraw = false, transparency = true, color = (:black, 0.1))
        streamplot!(ax1, field, -1..1, -1..1, colormap = cmap, gridsize= (40,40),
            arrow_size = 0.05,linewidth=1, transformation = (:xy, -zmax))
        streamplot!(ax2, field, -1..1, -1..1, colormap = cmap, gridsize= (40,40),
                arrow_size = 0.03,linewidth=1)
        fig[1,1] = ax1
        fig[1,2] = ax2
        fig
    end
    fig = with_theme(plot, theme_dark())
end

## 
function cubePlot()
    limits = Node(Rect(Vec3f0(-1), Vec3f0(2)))
    fig = Figure(resolution = (850, 800))
    ax = Axis3(fig, aspect = :data, perspectiveness = 0.5, xlabel = "time",
        ylabel = "longitude", zlabel = "latitude")
    cbox = cube!(ax; a1 = 1, a2 = 1,
        limits = limits, cmap = :CMRmap, transparent=false)
    ax.yticks = (-1:2:1,string.([lonRange[1], lonRange[2]]))
    ax.zticks = (-1:2:1,string.([latRange[1], latRange[2]]))
    ax.xticks = (-1:2:1,string.([timeRange[1],timeRange[2]]))
    cbar  = Colorbar(fig, cbox,  label = "temperature", height = Relative(0.5))
    #limits[] = Rect(Vec3f0(-2,-2,-2), Vec3f0(4.5))
    fig[1, 1] = ax
    fig[1, 2] = cbar
    fig
end



fig = Figure(resolution = (1200,600))
ax1 = Axis3(fig, aspect = (1,1,1), perspectiveness = 0.5,
            elevation = π/3.5, azimuth = 0.1π, )
volume(ax1, cube, algorithm = :absorption,
        transparency=true,colormap = cmap,
        figure=(backgroundcolor=:white) )













## =============================================================================
## ======== check the bistability region in 3d phase space =========
function check_bistability(param)
    output = [];
    # plt = plot()
    for i in 1:100
        u0 = rand(1:22., length(de1.states))
        prob = ODEProblem(de1, u0, tspan, param)
        sol = solve(prob)
        push!(output, sol[m1_HKCI][end]);
        # plot!(plt, sol, vars = [m1_HKCI])
    end 
    # display(plt)
    output = round.(output, digits=2)
    u_out= unique(output)
    out_min = minimum(output); out_max = maximum(output)
    separate_line = (out_max + out_min)/2
    left_bin = !isempty(output[output .< separate_line])
    right_bin = !isempty(output[output .> separate_line])
    bistable = left_bin && right_bin
    # bistable ? println("bistability region") : println("Non bistability region")
    # display(histogram(round.(output, digits=2), bins = 20, title = "bistability: $u_out"))
    return output, bistable
end 
output = check_bistability([1.5, 0.002, 0.991, 7.0, 0.025, 0.025, 0.0])

df_bistable_3d = DataFrame(K = Float64[], n = Float64[], up = Float64[], cri = [])
@showprogress for K = 0.011: 0.02: 1.0, n = 1.5:0.1:10.0, up = 1:10
    param = [up, 0.002, K, n, 0.025, 0.025, 0.0] # the last values if the output, make sure it is zero
    check_bistability(param)
    _, bistable = check_bistability(param)
    @show bistable
    record = [K, n, up, bistable]
    bistable ? push!(df_bistable_3d, record) : push!(df_bistable_3d, [K, n, up, false])
end 

sum(df_bistable_3d.cri)
CSV.write("./DEmodels/scripts_modi/3d_bistable.csv",df_bistable_3d)


df_bistable_3d  = CSV.File("./DEmodels/scripts_modi/3d_bistable.csv") |> DataFrame
plt_bistable_3d = @df df_bistable_3d scatter(
                    :n,
                    :K,
                    :up,
                    xlabel = "n",
                    ylabel = "K",
                    zlabel = "up",
                    label = ["Non bistable" "bistable"],
                    markersize = 2.2,
                    group =:cri,
                    dpi = 400)


## four different types for bistability and oscillation
df_bistable_3d
df_bifur_3d
rename!(df_bistable_3d, :cri => :bistability)
function merge_bistability_bifur(df_bifur_3d, df_bistable_3d)
end

combined = leftjoin(df_bistable_3d,df_bifur_3d,  on = [:K, :n,:up])
gdf = groupby(combined, [:bistability, :cri])
transform(gdf, [:bistability, :cri] => ((p, s) -> (claasification = string(p)*"_"*string(s))))


combined[:,:region] .= string(combined.bistability)
combined
CSV.write("./DEmodels/scripts_modi/3d_bifur_bistable_combined.csv",combined)





combined  = CSV.File("./DEmodels/scripts_modi/3d_bifur_bistable_combined.csv") |> DataFrame
# plot for combined 
plt_combine_3d = @df combined scatter(
                    :n,
                    :K,
                    :up,
                    xlabel = "n",
                    ylabel = "K",
                    zlabel = "up",
                    label = ["not bistable & no oscillation" "bistable & no oscillation" "bistable & oscillation"],
                    markersize = 2.2,
                    markercolor = ["blue" "orange" "green"],
                    group =:region,
                    dpi = 400)

## ======== slice layer plots 
layer_plot_set = []
for up_value  =  1.0 :  4.0
    layer_up = filter(row -> row.up == up_value, combined)
    layer_plt =@df layer_up scatter(
                        :n,
                        :K,
                        xlabel = "n",
                        ylabel = "K",
                        label = ["not bistable & no oscillation" "bistable & no oscillation" "bistable & oscillation"],
                        markersize = 4.8,
                        title = "up = $up_value",
                        group =:region,
                        markercolor = ["blue" "orange" "green"],
                        dpi = 400)
    push!(layer_plot_set,layer_plt)
end 
for up_value  =  5.0 :  10.0
    layer_up = filter(row -> row.up == up_value, combined)
    layer_plt =@df layer_up scatter(
                        :n,
                        :K,
                        xlabel = "n",
                        ylabel = "K",
                        label = ["bistable & no oscillation" "bistable & oscillation"],
                        markersize = 4.8,
                        title = "up = $up_value",
                        group =:region,
                        markercolor = [ "orange" "green"],
                        dpi = 400)
    push!(layer_plot_set,layer_plt)
end 

plot(layer_plot_set[[1,2,5,6,9,10]]..., layout = (3,2))
plot!(size=(1050,1050))
savefig("./DEmodels/scripts_modi/Vis/3d_bifur_bistable_combined_layer.png")  