using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

@variables x(t)
@parameters c(t) = 1

ev = ModelingToolkit.SymbolicDiscreteCallback(
    1.0 => [c ~ 1], discrete_parameters = c, iv = t)
@mtkcompile sys = System(
    D(x) ~ c * cos(x), t, [x], [c]; discrete_events = [ev])

prob = ODEProblem(sys, [x => 0.0], (0.0, 2pi))
sol = solve(prob, Tsit5())
sol[c]

bcl=1.0
duty=0.1
tstart=0
tend=10
ts0 = collect(tstart:bcl:tend)
ts1 = ts0 .+ duty
ev0 = ModelingToolkit.SymbolicDiscreteCallback(
    ts0 => [c ~ 3], discrete_parameters = c, iv = t)
ev1 = ModelingToolkit.SymbolicDiscreteCallback(
    ts1 => [c ~ 1], discrete_parameters = c, iv = t)

@mtkcompile sys = System(D(x) ~ c * cos(x), t; discrete_events = [ev0, ev1])

discrete_events(sys)
