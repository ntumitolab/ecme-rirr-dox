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
@named sys = build_model(; tend, bcl)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
prob = ODEProblem(sys, [u0; sys.kf16_C1 => 2Hz / μM; sys.kf17_C1 => 0.02Hz / μM; sys.K010_C3 => 3Hz / mM; sys.KCAT_IDH => 43Hz], tend)

@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true)

for i in sts
    istr = replace(string(i), "(t)" => "")
    println("sys.", istr, " => ", sol[i][end], ",")
end

#---
plot(sol, idxs=sys.vm, legend=:right, tspan=(900second, 901second)) |> PNG

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
plot(sol, idxs = [sys.sox_i, sys.sox_m], tspan=(900e3, 910e3))

# ROS
plot(sol, idxs = [sys.vROSIf, sys.vROSIq, sys.vROSC1, sys.vROSC3])

# O2 Shunt
plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="O2 Shunt")

#---
plot(sol, idxs = [sys.dpsi], tspan=(900e3, 910e3))

#---
plot(sol, idxs = [sys.vC5], tspan=(900e3, 910e3))
