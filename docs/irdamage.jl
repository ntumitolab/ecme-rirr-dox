# Ischemic reperfusion damage
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using DiffEqCallbacks
using ECMEDox
using ECMEDox: second, Hz, μM, nM
using Plots
using DisplayAs: PNG

#---
tend = 400.0second
bcl = 1.0second
@named sys = build_model(; bcl, tend)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
@unpack O2 = sys
prob = ODEProblem(sys, u0, tend, [O2 => 6nM])

affect!(integrator) = integrator.ps[O2] = 6μM
cb = PresetTimeCallback(60.0second, affect!)

#---
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true, dt=1e-6, callback=cb)

@unpack vm, dpsi, atp_i, adp_i = sys
plot(sol, idxs=dpsi)

plot(sol, idxs=vm)

plot(sol, idxs=[atp_i/ adp_i]) |> PNG

@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal= sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], legend=:right, title="CAC metabolites") |> PNG

# Because Qdot_n recovers?
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol, idxs=[Q_n + Q_p, Qdot_n, QH2_n+QH2_p, Qdot_p], title="Q cycle", legend=:right) |> PNG

plot(sol, idxs=[Q_n + Q_p, Qdot_n, QH2_n+QH2_p, Qdot_p], title="Q cycle", legend=:right, tspan=(185second, 190second)) |> PNG

plot(sol, idxs=[sys.vHresC1, sys.vHresC3, sys.vHresC4], title="ETC resp rate", legend=:right) |> PNG

# Q cycle revive here
plot(sol, idxs=[Q_n, Q_p], title="Q", legend=:right, tspan=(180second, 184second)) |> PNG

plot(sol, idxs=[fes_ox, fes_rd, cytc_ox, cytc_rd], title="Q cycle", legend=:right) |> PNG

plot(sol, idxs=[sys.vROSC1, sys.vROSC3], legend=:right, title="Superoxide generation") |> PNG

plot(sol, idxs=[sys.sox_i, sys.sox_m], legend=:right, title="Superoxide") |> PNG

plot(sol, idxs=[sys.vQC1], legend=:right, title="Complex I rates") |> PNG

@unpack vCS, vACO, vIDH, vKGDH, vSL, vFH, vMDH, vAAT, vSDH = sys
plot(sol, idxs=[vCS, vACO, vIDH, vKGDH, vFH, vMDH, vAAT, vSDH], legend=:right, title="CAC flux") |> PNG

# O2 shunt
idxs = (sys.t/1000, sys.vROS / (sys.vO2 + sys.vROS))
fig = plot(sol, idxs=idxs, title="ROS shunt fraction") |> PNG
