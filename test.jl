using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DiffEqCallbacks
using Plots

@parameters begin
    kf = 1
    keq = 1
    kb = kf / keq
end

@variables a(t) = 1 b(t) = 0

v = kf * a - kb * b

@mtkcompile sys = System([D(a) ~ -v; D(b) ~ v], t)

prob = ODEProblem(sys, [], (0.0, 1.0))
