#===
# Doxorubicin addition test
===#
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DataFrames
using CSV
using ECMEDox
using ECMEDox: second, mM, Hz, μM
Plots.default(lw=2, size=(600, 600))

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX, ρC4, ρC3 = sys
sts = unknowns(sys)
obs = [i.lhs for i in observed(sys)]
u0 = build_u0(sys)

# The model is very sensitive to aconitase activity
# The phase transition (of the Q cycle) is between 250uM to 260uM of DOX
prob = ODEProblem(sys, u0, tend)
alg = FBDF()
@time sol = solve(prob, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

#---
plot(sol, idxs=sys.vm, tspan=(100second, 105second), plotdensity=500)
#---
plot(sol, idxs=[sys.atp_m / sys.adp_m], plotdensity=1000)
#---
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal])
# Calcium transient too low, weird ICa flux?
plot(sol, idxs=[sys.ca_i, sys.ca_m], tspan=(100second, 105second))
plot(sol, idxs=[sys.ICaL, sys.INaCa], tspan=(100second, 105second))

@unpack c0_lcc, c1_lcc, c2_lcc, c3_lcc, c4_lcc, o_lcc, cca0_lcc, cca1_lcc, cca2_lcc, cca3_lcc, cca4_lcc, x_yca = sys
plot(sol, idxs=[o_lcc, x_yca], tspan=(100second, 102second))
plot(sol, idxs=[sys.ICaK], tspan=(100second, 102second))
plot(sol, idxs=[sys.vUni, sys.vNaCa], tspan=(100second, 105second))

prob0 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC4 => 325μM])
prob1 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC4 => 500μM])
prob2 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC3 => 500μM])
alg = FBDF()
@time sol0 = solve(prob0, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)
@time sol1 = solve(prob1, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)
@time sol2 = solve(prob2, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

@unpack atp_i, adp_i, vm, na_o, na_i, atp_i, adp_i, atp_m, adp_m, sox_i, sox_m = sys
@unpack ca_i, ca_nsr, ca_jsr, ca_ss = sys
plot(sol0, idxs=vm, label="C4 1x", title="PM potential")
plot(sol0, idxs=[atp_i / adp_i], label="C4 1x", title="ATP:ADP", density=1000)

plot(sol1, idxs=vm, label="C4 3x")
plot(sol2, idxs=vm, label="C3 3x")

p1 = plot(sol, idxs=[ca_i, ca_ss], tspan=(990, tend))
p2 = plot(sol, idxs=[ca_nsr, ca_jsr], tspan=(990, tend))
plot(p1, p2, layout=(2, 1))
@unpack IK1, IK, IKp, IKatp, INa, INaCa, ICaL = sys
plot(sol0, idxs=[IK1, IK, IKp, IKatp, ICaL], tspan=(400second, tend))

plot(sol0, idxs=[atp_i / adp_i], label="C4 1x", title="ATP:ADP")
plot!(sol1, idxs=[atp_i / adp_i], label="C4 3x", line=:dash)
plot!(sol2, idxs=[atp_i / adp_i], label="C3 3x", line=:dot)

@unpack vHres, vHu, vANT, vHleak, vNaCa, vUni, vIMAC, vROS, vTrROS, vO2 = sys
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol0, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="C4 1x")
plot(sol1, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="C4 3x")
plot(sol2, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="C3 3x")

plot(sol0, idxs=[sys.dpsi * 1000], title="C4 1x", tspan=(700, 800))
plot!(sol1, idxs=[sys.dpsi * 1000], title="C4 3x", tspan=(700, 800))
plot!(sol2, idxs=[sys.dpsi * 1000], title="C3 3x", ylims=(100, 160), tspan=(700, 800))

plot(sol, idxs=[vHres, vHu, vANT, vHleak, vNaCa, vUni, vTrROS], ylims=(0, 5))
plot(sol, idxs=vANT, tspan=(400, 500))
plot(sol, idxs=[fes_ox / fes_rd], tspan=(400, 500))
plot(sol, idxs=[cytc_ox / cytc_rd], tspan=(400, 500))
plot(sol, idxs=[sys.sox_i, sys.h2o2_i], tspan=(740, 800))
plot(sol, idxs=[sys.vSOD_i, sys.vGPX_i, sys.vGR_i, sys.vCAT], tspan=(400, 500))
plot(sol0, idxs=[sys.vROSC1, sys.vROSC3], tspan=(700, 800), title="C4 1x")
plot(sol1, idxs=[sys.vROSC1, sys.vROSC3], tspan=(700, 800), title="C4 3x")
plot(sol2, idxs=[sys.vROSC1, sys.vROSC3], tspan=(700, 800), title="C3 3x")

plot(sol0, idxs=vO2, label="C4 1x", title="vO2")
plot!(sol1, idxs=vO2, label="C4 3x", line=:dash)
plot!(sol2, idxs=vO2, label="C3 3x", line=:dot)

plot(sol0, idxs=vROS, label="C4 1x", title="vROS")
plot!(sol1, idxs=vROS, label="C4 3x")

plot(sol, idxs=[sys.vROS / (sys.vO2 + sys.vROS)])

@unpack nad_m, nadh_m, vC1, vIDH, vKGDH, vMDH = sys
plot(sol, idxs=[nad_m / nadh_m])

@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
@unpack vCS, vACO, vIDH, vKGDH, vSL, vFH, vMDH, vAAT, vSDH = sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal])
