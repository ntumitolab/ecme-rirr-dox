# # Initial conditions
# 1Hz pacing for 1000 seconds.
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ECMEDox
using ECMEDox: second, mM, Hz, μM, mV
using Plots
#---
tend = 1000.0second
bcl = 1.0second
@time @named sys = build_model()
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()

for eq in equations(sys)
    println(eq)
end

@unpack iStim = sys
callback = build_stim_callbacks(iStim, tend; period=bcl)
prob = ODEProblem(sys, u0, tend)

#---
@time sol = solve(prob, alg; callback = callback, saveat=0.01second);
