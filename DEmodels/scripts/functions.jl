using DifferentialEquations


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
