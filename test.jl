using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DiffEqCallbacks
using Plots

@parameters iStim(t)
@variables x(t) v_istim(t)

@named sys = ODESystem([D(x) ~ -x + iStim], t)
sys = structural_simplify(sys)

stimulate! = (integrator) -> begin
    integrator.ps[iStim] = 5
end

tspan = (0.0, 10.0)
prob = ODEProblem(sys, [x => 0.0, iStim=>0.0], tspan)
cb = PresetTimeCallback(1.0, stimulate!)

@time sol = solve(prob; callback=cb)
