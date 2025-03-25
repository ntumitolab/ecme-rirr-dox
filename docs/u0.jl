# # Initial conditions
# 1Hz pacing for 1000 seconds
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using ECMEDox
using ECMEDox: second, mM, Hz, μM, cm, μAcm⁻²
using Plots
using DisplayAs: PNG

tend = 1000.0second
bcl = 1.0second
@named sys = build_model(; bcl, tend)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
prob = ODEProblem(sys, u0, tend, [sys.k_IMAC => -0.07])
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true)

for i in sts
    istr = replace(string(i), "(t)" => "")
    println("sys.", istr, " => ", sol[i][end], ",")
end
