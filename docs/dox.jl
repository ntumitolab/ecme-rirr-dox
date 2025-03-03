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
using ECMEDox: second, μM

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX, ρC4, ρC3 = sys
sts = unknowns(sys)
u0 = build_u0(sys)

# ## No DOX (control)
prob = ODEProblem(sys, u0, tend, [DOX => 0μM])
alg = Rodas5P()
opts = (; reltol=1e-6, abstol=1e-6, progress=true, maxiters=1e8)
@time sol = solve(prob, alg; opts...)

# The phase transition (of the Q cycle) is between 250uM to 260uM of DOX
prob = ODEProblem(sys, u0, tend, [DOX => 250μM])
alg = FBDF()
opts = (; reltol=1e-6, abstol=1e-6, progress=true, maxiters=1e8)
@time sol = solve(prob, alg; opts...)
idxs = (sys.t/1000, sys.vm)
plot(sol, idxs=idxs, tspan=(100second, 103second), lab=false, title="PM potential")

#---
idxs = (sys.t/1000, [sys.atp_i / sys.adp_i])
plot(sol, idxs=idxs, tspan=(100second, 103second), lab="DOX=250", title="ATP:ADP")

#---
prob0 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC4 => 325μM])
prob1 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC4 => 500μM])
prob2 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC3 => 500μM])
@time sol0 = solve(prob0, alg; opts...)
@time sol1 = solve(prob1, alg; opts...)
@time sol2 = solve(prob2, alg; opts...)

idxs = (sys.t/1000, sys.vm)

#---
plot(sol0, idxs=idxs, lab="DOX=260", title="PM potential")

#---
plot(sol1, idxs=idxs, lab="DOX=260, ρC4=500", title="PM potential")

#---
plot(sol2, idxs=idxs, lab="DOX=260, ρC3=500", title="PM potential")

#---
idxs = (sys.t/1000, [sys.atp_i / sys.adp_i])
tspan=(100second, 200second)
plot(sol0, idxs=idxs, label="C4 325", title="ATP:ADP"; tspan)
plot!(sol1, idxs=idxs, label="C4 500", line=:dash; tspan)
plot!(sol2, idxs=idxs, label="C3 500", line=:dot; tspan)

#---
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol0, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="Q cycle", legend=:right)

#---
idxs = (sys.t/1000, sys.dpsi)
tspan=(100second, 200second)
plot(sol0, idxs=idxs, title="MMP", lab="C4 325" ; tspan)
plot!(sol1, idxs=idxs, label="C4 500"; tspan)
plot!(sol2, idxs=idxs, label="C3 500"; tspan, xlabel="Time (s)", ylabel="MMP (mV)")

#---
idxs = (sys.t/1000, sys.vO2)
plot(sol0, idxs=idxs, label="C4 325", title="Oxygen consumption")
plot!(sol1, idxs=idxs, label="C4 500", line=:dash)
plot!(sol2, idxs=idxs, label="C3 500", line=:dot, xlabel="Time (s)", ylabel="vO2 (μM/ms)")

#---
idxs = (sys.t/1000, sys.vROS / (sys.vO2 + sys.vROS))
plot(sol0, idxs=idxs, label="C4 325", title="ROS vO2 fraction")
plot!(sol1, idxs=idxs, label="C4 500", line=:dash)
plot!(sol2, idxs=idxs, label="C3 500", line=:dot, xlabel="Time (s)")

#---
idxs = (sys.t/1000, sys.nadh_m)
plot(sol0, idxs=idxs, label="C4 325", title="NADH")
plot!(sol1, idxs=idxs, label="C4 500", line=:dash)
plot!(sol2, idxs=idxs, label="C3 500", line=:dot, xlabel="Time (s)")
