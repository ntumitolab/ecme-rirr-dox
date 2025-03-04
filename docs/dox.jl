#===
# Doxorubicin addition test
===#
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DisplayAs
using ECMEDox
using ECMEDox: second, mM, Hz, μM
Plots.default(lw=2, size=(600, 600))

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX, ρC4, ρC3, KEQ_C5 = sys
sts = unknowns(sys)
u0 = build_u0(sys)

# The phase transition (of the Q cycle) is between 240uM to 260uM of DOX
prob = ODEProblem(sys, u0, tend, [DOX => 240μM])
alg = FBDF()
opts = (; reltol=1e-6, abstol=1e-6, progress=true, maxiters=1e8)
@time sol = solve(prob, alg; opts...)
idxs = (sys.t/1000, sys.vm)
plot(sol, idxs=idxs, lab=false, title="PM potential")
#---
idxs = (sys.t/1000, sys.dpsi)
plot(sol, idxs=idxs, lab=false, title="Mito potential")
#---
idxs = (sys.t/1000, [sys.atp_i / sys.adp_i])
plot(sol, idxs=idxs, lab="DOX=240", title="ATP:ADP")

#---
prob0 = ODEProblem(sys, u0, tend, [DOX => 250μM])
prob1 = ODEProblem(sys, u0, tend, [DOX => 250μM, ρC4 => 500μM])
prob2 = ODEProblem(sys, u0, tend, [DOX => 250μM, ρC3 => 500μM])
@time sol0 = solve(prob0, alg; opts...)
@time sol1 = solve(prob1, alg; opts...)
@time sol2 = solve(prob2, alg; opts...)

idxs = (sys.t/1000, sys.vm)

#---
plot(sol0, idxs=idxs, lab="DOX=250", title="PM potential")

#---
plot(sol1, idxs=idxs, lab="DOX=250, ρC4=500", title="PM potential")

#---
plot(sol2, idxs=idxs, lab="DOX=250, ρC3=500", title="PM potential")

#---
idxs = (sys.t/1000, [sys.atp_i / sys.adp_i])
plot(sol0, idxs=idxs, label="DOX=250", title="ATP:ADP")
plot!(sol1, idxs=idxs, label="C4 500uM")
plot!(sol2, idxs=idxs, label="C3 500uM")

#---
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol0, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="Q cycle", legend=:right)

#---
idxs = (sys.t/1000, sys.vO2 + sys.vROS)
plot(sol0, idxs=idxs, label="DOX=250", title="Oxygen consumption")
plot!(sol1, idxs=idxs, label="C4 500")
plot!(sol2, idxs=idxs, label="C3 500", xlabel="Time (s)", ylabel="vO2 (μM/ms)")

#---
idxs = (sys.t/1000, sys.vROS / (sys.vO2 + sys.vROS))
plot(sol0, idxs=idxs, label="DOX=250", title="ROS vO2 shunt fraction")
plot!(sol1, idxs=idxs, label="C4 500")
plot!(sol2, idxs=idxs, label="C3 500", xlabel="Time (s)")

# Succinate accumulation adn CAC slowed down.
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal= sys
plot(sol0, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], leg=:right, title="CAC metabolites")

# NADH production from CAC decreased
idxs = (sys.t/1000, sys.nadh_m)
plot(sol0, idxs=idxs, label="DOX=250", title="NADH (mito)")
plot!(sol1, idxs=idxs, label="C4 500")
plot!(sol2, idxs=idxs, label="C3 500", xlabel="Time (s)")
