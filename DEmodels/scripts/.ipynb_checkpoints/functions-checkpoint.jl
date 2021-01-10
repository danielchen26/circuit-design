# using DifferentialEquations
using OrdinaryDiffEq
using Optim



# =================== Using this right now =================

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




# ======= control problems initialization ===========
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


#  ========================= Plotting functions =================================
function localminmax_sig(sol1, t_on, opt)
    plt6_min = plot(sol1, vars=[:m1_HKCI], color =:green, xlims = (1:t_on))
    locs6_min =  findlocalminima(sol1[6,:])
    mark6_min = [i[1] for i in locs6_min]
    scatter!(sol1[mark6_min], vars=[:m1_HKCI], xlims = (1:t_on), marker = (:circle, 1, 0.6, :purple, stroke(2, 0.5, :orange, :dot)))
    title!("Gate 6")

    plt6_max = plot(sol1, vars=[:m1_HKCI], color =:green, xlims = (1:t_on))
    locs6_max =  findlocalmaxima(sol1[6,:])
    mark6_max = [i[1] for i in locs6_max]
    scatter!(sol1[mark6_max], vars=[:m1_HKCI], xlims = (1:t_on), marker = (:circle, 1, 0.6, :purple, stroke(2, 0.5, :orange, :dot)))
    title!("Gate 6")



    plt7_min = plot(sol1, vars=[:m1_PhlF], color = :darkorange, xlims = (1:t_on))
    locs7_min =  findlocalminima(sol1[7,:])
    mark7_min = [i[1] for i in locs7_min]
    scatter!(sol1[mark7_min], vars=[:m1_PhlF], xlims = (1:t_on), marker = (:circle, 1, 0.6, :red, stroke(2, 0.5, :blue, :dot)))
    title!("Gate 7")

    plt7_max = plot(sol1, vars=[:m1_PhlF], color = :darkorange, xlims = (1:t_on))
    locs7_max =  findlocalmaxima(sol1[7,:])
    mark7_max = [i[1] for i in locs7_max]
    scatter!(sol1[mark7_max], vars=[:m1_PhlF], xlims = (1:t_on), marker = (:circle, 1, 0.6, :red, stroke(2, 0.5, :blue, :dot)))
    title!("Gate 7")

    plt_min = plot(plt6_min,plt7_min, layout = (2,1))
    plt_max = plot(plt6_max,plt7_max, layout = (2,1))

    if opt == "lmin"
        return plt_min, mark6_min, mark7_min
    elseif opt == "lmax"
        return plt_max, mark6_max, mark7_max
    end
end



# Make animation for control problems ======
function Anim_gen_controlproblem(cycle, Δ0,  Δ,  δ, A, sig_rg)
    anim = @animate for δ in sig_rg
        time, signal = signal_gen(cycle, Δ0,  Δ,  δ, A)
        ts, cb = cb_gen([time...], signal...)
        sol = solve(prob0, Tsit5(), callback=cb, tstops=ts, reltol = 1e-15, abstol = 1e-19)
        py = plot(sol,vars=[:m1_HKCI,:m1_PhlF], lw =2,xlabel = "time", ylabel = "concentration", title = "Signal Duration: $δ")
    end
    return anim
end



# === Check whether the gate switches or not ============
function switching(sol, cycles)
    # println("End of relaxation time: ", Δ0)
    G6 = sol(Δ0)[6];  G7 = sol(Δ0)[7]
    # @show G6, G7

    g6ii_tot = g7ii_tot = 0
    g6i_tot = g7i_tot = 0
    if G6 > G7
        # locs = findlocalmaxima(sol[7, :])
        # marklist = [i[1] for i in locs]
        # marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
        # display(scatter!(sol[marklist], vars = [:m_PhlF], xlims = (1:2000.0), marker = marker_spec)) # xlims = (1:t_on),
        # require cycle(i+1)_mid/cycle(i)_mid = up/dn for both gate 6 & 7
        # t0 = [Δ0, Δ0 + δ]
        # i = 1 # default
        # tii = Δ0 + i * Δ
        # ti = Δ0 + (i - 1) * Δ
        #
        # @show tii, ti
        # @show g6ii = sol(tii)[6]
        # @show g6i = sol(ti)[6]
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

            # R_g6 = g6ii/g6i
            # R_g7 = g7ii/g7i
            # R_accum *= R_g7*R_g6
            # # switch_sign = sign((g6ii -g6i)*(g7ii -g7i))
            # @show R_g6, R_g7, R_accum
        end
    elseif G6 < G7
        # locs = findlocalmaxima(sol[6, :])
        # marklist = [i[1] for i in locs]
        # marker_spec = (:circle, 5, 0.6, :purple, stroke(6, 0.5, :orange, :dot))
        # display(scatter!(sol[marklist], vars = [:m_HKCI], xlims = (1:2000.0), marker = marker_spec)) # xlims = (1:t_on),
        # require cycle(i+1)_mid/cycle(i)_mid = up/dn for both gate 6 & 7
        # t0 = [Δ0, Δ0 + δ]
        # i = 1 # default
        # tii = Δ0 + i * Δ
        # ti = Δ0 + (i - 1) * Δ
        #
        # @show tii, ti
        # @show g6ii = sol(tii)[6]
        # @show g6i = sol(ti)[6]

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








# # =============== deprecated function =============================
# function cb_gen(ts, vars...)
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#         for i in eachindex(ts)
#             if integrator.t == ts[i]
#                 integrator.p = vars[i]
#             end
#         end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     # @show vars
#     return ts, cb
# end

# # ======= Callback functions ==========
# function make_cb(ts_in, vars...)
#     ts = ts_in
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#         for  i = 1:2: length(vars)
#             if integrator.t == ts[1]
#                 integrator.p[vars[i]] = vars[i+1]
#             elseif integrator.t == ts[2]
#                 integrator.p[vars[i]] = 0.0
#             end
#         end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     @show vars
#     return ts, cb
# end
#
# function make_cb2(ts, vars...)
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#       if integrator.t == ts[1]
#           integrator.p[1] = vars[1]
#       elseif integrator.t == ts[2]
#           integrator.p[1] = vars[2]
#       elseif integrator.t == ts[3]
#           integrator.p[1] = vars[3]
#       elseif integrator.t == ts[4]
#           integrator.p[1] = vars[4]
#       end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     @show vars
#     return ts, cb
# end
# # make_cb2([1000,1500,3000,3600],  20., 0., 20., 0.)
#
# function mk_cb_mul(ts, vars...)
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#         for i in eachindex(ts)
#             if integrator.t == ts[i]
#                 integrator.p[1] = vars[i]
#             end
#         end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     @show vars
#     return ts, cb
# end
# mk_cb_mul([300,400,500,600], 0., 20., 0., 20.)
