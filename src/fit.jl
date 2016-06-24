include("mueller.jl")
using NLopt

function cost_circ(x, grad, data)
    stokes = Array{Float64}(size(data,1), 4)
    for (i, d) in enumerate(data)
        stokes[i,:] = müller_rec()*
        müller_sky(d.pa_az)*
        s_source(flux = 0.0, pa = 0.0)
        
# Set up NLopt
opt = Opt(:LN_COBYLA, nfields(FlareModel) - 3)
xtol_rel!(opt, σ/10)

initial_step!(opt, init_steps)
lower_bounds!(opt, lower_bounds)
upper_bounds!(opt, upper_bounds)

min_objective!(opt, (x, grad) -> cost(x, grad, data, vec(weights[j,:]), FlareModel())) # model, point
minf, nlopt_return = optimize!(opt, flare)

