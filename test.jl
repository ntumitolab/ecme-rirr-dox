using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DiffEqCallbacks
using Plots

@parameters iStim(t)
@variables x(t) v_istim(t)

@named sys = ODESystem([D(x) ~ -x + v_istim, v_istim ~ iStim], t)
simp = structural_simplify(sys)

setter! = ModelingToolkit.setp(simp, iStim)
stimulate! = (integrator) -> begin
    setter!(integrator, 5)
end

tspan = (0.0, 10.0)
prob = ODEProblem(simp, [x => 0.0], tspan, [iStim=>0.0])
cb = PresetTimeCallback(1.0, stimulate!)

@time sol = solve(prob; callback=cb)

sol[x]
sol[v_istim]

plot(sol, idxs=[x, v_istim])

parameter_values(prob)
