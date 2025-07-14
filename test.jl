using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
using DiffEqCallbacks

@variables a(t) = 1 b(t) = 0
@parameters kf(t) = 1 kb(t) keq(t)=1

v = kf * a - kb * b
eqs = [
    D(a) ~ -v,
    D(b) ~ v,
    kb ~ kf / keq
]

@mtkcompile sys = System(eqs, t)

prob = ODEProblem(sys, [], (0.0, 10.0))

prob.ps

sol = solve(prob)
sol.ps[kb]
plot(sol, idxs=[a, b, kb])

prob2 = remake(prob)
prob2.ps[keq] = 10.0
prob2.ps[kb]

sol2 = solve(prob2)
plot(sol2, idxs=[a, b, kb])

event! = (integrator) -> begin
    integrator.ps[sys.keq] = 10
end

callback = PresetTimeCallback(5.0, event!)

sol3 = solve(prob; callback)

plot(sol3, idxs=[a, b, kb])
