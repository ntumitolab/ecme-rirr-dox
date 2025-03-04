# # Initial conditions
# 1Hz pacing for 1000 seconds
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using ECMEDox
using ECMEDox: second, mM, Hz, μM

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
sts = unknowns(sys)
u0 = build_u0(sys)

# The model is very sensitive to aconitase activity
# The phase transition (of the Q cycle) is between 250uM to 260uM of DOX
prob = ODEProblem(sys, u0, tend)
alg = FBDF()
@time sol = solve(prob, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8, save_everystep=false)

for i in sts
    println("sys.", i, " => ", sol[i][end], ",")
end
