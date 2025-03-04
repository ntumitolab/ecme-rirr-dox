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

prob = ODEProblem(sys, u0, tend, [sys.P_CA_LCC => 8e-4cm * Hz])
alg = KenCarp4()
@time sol = solve(prob, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

for i in sts
    println("sys.", i, " => ", sol[i][end], ",")
end
