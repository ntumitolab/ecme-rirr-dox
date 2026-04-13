# # Initial conditions
# 1Hz pacing for 1000 seconds.
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Model
using Model: second, mM, Hz, μM, mV
using Plots
#---
tend = 1000.0second
bcl = 1.0second
@time @named sys = build_model()
u0 = build_u0(sys)
sts = unknowns(sys)
alg = FBDF()

for eq in equations(sys)
    println(eq)
end

@unpack iStim = sys
callback = build_stim_callbacks(iStim, tend; period=bcl)
prob = ODEProblem(sys, u0, tend)

#---
@time sol = solve(prob, alg; callback = callback, saveat=0.01second);

for i in sts
    istr = replace(string(i), "(t)" => "")
    println("sys.", istr, " => ", sol[i][end], ",")
end

# ## Action potential
plot(sol, idxs=sys.vm, legend=:right, tspan=(900second, 901second))

# ## Citric acid cycle metabolites
# Citrate and isocitrate have the highest concentrations.
@unpack citrate, isocitrate, oaa, akg, scoa, succinate, fumarate, malate = sys
plot(sol, idxs=[citrate, isocitrate, oaa, akg, scoa, succinate, fumarate, malate], legend=:right, title="CAC metabolites")

# ## CAC flux
@unpack vMDH, vAAT, vIDH = sys
plot(sol, idxs=[vMDH, vAAT, vIDH], legend=:right, title="CAC flux")

# ## Q cycle
@unpack Q_n, SQn, QH2_n, QH2_p, Q_p, SQp, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol, idxs=[Q_n + Q_p, SQn, QH2_n + QH2_p, SQp], title="Q cycle", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")

#---
plot(sol, idxs=[fes_ox, fes_rd, cytc_ox, cytc_rd], title="Q cycle (downstream)", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")

# ## Proton pumping
plot(sol, idxs = [sys.vHresC1, sys.vHresC3, sys.vHresC4], ylims=(0, 3))

# ## ROS
plot(sol, idxs = [sys.sox_i, sys.sox_m], tspan=(900e3, 910e3))

#---
plot(sol, idxs = [sys.vROSIf, sys.vROSIq, sys.vROSC1, sys.vROSC3], tspan=(900e3, 910e3))

#---
plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="O2 Shunt", tspan=(900e3, 910e3), ylims=(0, 5))

# ## MMP

plot(sol, idxs = [sys.dpsi], tspan=(900e3, 910e3))

#---
plot(sol, idxs = [sys.vC5], tspan=(900e3, 910e3))
