# # Initial conditions
# 1Hz pacing for 1000 seconds
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using ECMEDox
using ECMEDox: second, mM, Hz, μM, build_stim_callbacks
using Plots
using DisplayAs: PNG

tend = 1000.0second
bcl = 1.0second
@named sys = build_model(; tend, bcl=0)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
prob = ODEProblem(sys, u0, tend)

cb = build_stim_callbacks(sys.iStim, tend)
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true, callback=cb)

for i in sts
    istr = replace(string(i), "(t)" => "")
    println("sys.", istr, " => ", sol[i][end], ",")
end

# Citric acid cycle metabolites
# Dominated by citrate and isocitrate
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], legend=:right, title="CAC metabolites") |> PNG

# Q cycle
@unpack Q_n, SQn, QH2_n, QH2_p, Q_p, SQp, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
pl_q = plot(sol, idxs=[Q_n + Q_p, SQn, QH2_n + QH2_p, SQp], title="Q cycle", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)") |> PNG

#---
pl_q = plot(sol, idxs=[fes_ox, fes_rd, cytc_ox, cytc_rd], title="Q cycle (downstream)", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)") |> PNG

#---
plot(sol, idxs = [sys.vHresC1, sys.vHresC3, sys.vHresC4], ylims=(0, 3))

#---
plot(sol, idxs = [sys.sox_i, sys.sox_m])

# ROS
plot(sol, idxs = [sys.vROSC1, sys.vROSC3])

# O2 Shunt
plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="O2 Shunt")

#---
plot(sol, idxs = [sys.dpsi])

#---
plot(sol, idxs=sys.sox_i)

sol[sys.ΔVROS]
