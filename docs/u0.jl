# # Initial conditions
# 1Hz pacing for 1000 seconds
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using ECMEDox
using ECMEDox: second, mM, Hz, μM, cm
using Plots
using DisplayAs: PNG

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
sts = unknowns(sys)
u0 = build_u0(sys)

prob = ODEProblem(sys, u0, tend, [])
alg = KenCarp4()
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true)

for i in sts
    println("sys.", i, " => ", sol[i][end], ",")
end
