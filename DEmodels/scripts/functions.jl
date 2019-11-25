using DifferentialEquations

# ======= Callback functions ==========
function make_cb(ts_in, vars...)
    ts = ts_in
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
        for  i = 1:2: length(vars)
            if integrator.t == ts[1]
                integrator.p[vars[i]] = vars[i+1]
            elseif integrator.t == ts[2]
                integrator.p[vars[i]] = 0.0
            end
        end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    @show vars
    return ts, cb
end

function make_cb2(ts, vars...)
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
      if integrator.t == ts[1]
          integrator.p[1] = vars[1]
      elseif integrator.t == ts[2]
          integrator.p[1] = vars[2]
      elseif integrator.t == ts[3]
          integrator.p[1] = vars[3]
      elseif integrator.t == ts[4]
          integrator.p[1] = vars[4]
      end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    @show vars
    return ts, cb
end
# make_cb2([1000,1500,3000,3600],  20., 0., 20., 0.)

function mk_cb_mul(ts, vars...)
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
        for i in eachindex(ts)
            if integrator.t == ts[i]
                integrator.p[1] = vars[i]
            end
        end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    @show vars
    return ts, cb
end
# mk_cb_mul([300,400,500,600], 0., 20., 0., 20.)





# # =============== Multi-bit counter Cost function =============================
# using Parameters
#
# function cost_Mbit(P::Hill, Δ, δ, sol)
#     # First estimate the first oscillation time
#
#     # Test two adjacent oscillations
#     G6 = sol[6,end]; G7 =sol[7,end]
#
# end
