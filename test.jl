using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DisplayAs: PNG
using ECMEDox
using ECMEDox: second, mM, Hz, μM
Plots.default(lw=1.5, size=(600, 600))

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX, ET_SOD_I = sys
sts = unknowns(sys)
u0 = build_u0(sys)
alg = KenCarp47()
opts = (; reltol=1e-6, abstol=1e-6, progress=true)

prob = ODEProblem(sys, u0, tend)
prob310 = ODEProblem(sys, u0, tend, [DOX => 310μM, ET_SOD_I=>1μM])


sol = solve(prob, alg; opts...)  ## Warm-up
sol310 = solve(prob310, alg; opts...)
sol310[sys.dpsi]

plot(sol310, idxs=sys.dpsi)
plot(sol310, idxs=[sys.sox_i, sys.sox_m])
plot(sol310, idxs=[sys.vHres, sys.vHu, sys.vIMAC])
plot(sol310, idxs=sys.gsh_i)-+-+    -+  -+  -
