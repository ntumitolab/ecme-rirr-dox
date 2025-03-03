# # Initial conditions
# 1Hz pacing for 1000 seconds
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using ECMEDox
using ECMEDox: second, mM

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
sts = unknowns(sys)
u0 = build_u0(sys)

prob = ODEProblem(sys, u0, tend, [sys.ρC1 => 5mM])
alg = FBDF()
@time sol = solve(prob, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

for i in sts
    println("sys.", i, " => ", sol[i][end], ",")
end
