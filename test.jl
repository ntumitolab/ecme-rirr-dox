using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DiffEqCallbacks
using Plots

@parameters iStim(t)
@variables x(t) v_istim(t) rhs(t)

stimulate! = (integrator) -> begin
    integrator.ps[iStim] = 5
end
cb = PresetTimeCallback(1.0, stimulate!)

ev = ModelingToolkit.SymbolicDiscreteCallback(
    [1.0] => [iStim ~ 5], discrete_parameters = iStim, iv = t)

@mtkcompile sys = ODESystem([D(x) ~ rhs, rhs ~ -x + v_istim, v_istim ~ iStim], t; discrete_events = [ev])

tspan = (0.0, 10.0)
prob = ODEProblem(sys, [x => 0.0, iStim=>0.0], tspan)

@time sol = solve(prob)

ts = 0:0.1:10.0
ys = sol(ts, idxs=[sys.iStim, sys.x]).u |> stack |> transpose
plot(ts, ys, label=["iStim" "x"])
