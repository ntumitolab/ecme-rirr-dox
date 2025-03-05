#===
# Doxorubicin addition test
===#
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DisplayAs: PNG
using ECMEDox
using ECMEDox: second, mM, Hz, μM
Plots.default(lw=1.5, size=(600, 600))

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX, ρC4, ρC3 = sys
sts = unknowns(sys)
u0 = build_u0(sys)

# The phase transition (of the Q cycle) is between 293uM and 294uM of DOX
prob = ODEProblem(sys, u0, tend, [DOX => 295μM])
alg = KenCarp4()
opts = (; reltol=1e-6, abstol=1e-6, progress=true, maxiters=1e8)
@time sol = solve(prob, alg; opts...)

#---
s = sys.t/1000
@unpack vm = sys
plot(sol, idxs=(s, vm), lab=false, title="PM potential") |> PNG
#---
@unpack dpsi = sys
plot(sol, idxs=(s, dpsi), lab=false, title="Mito potential") |> PNG
#---
@unpack atp_i, adp_i = sys
plot(sol, idxs=(s, [atp_i / adp_i]), lab="DOX=295", title="ATP:ADP") |> PNG

#---
prob0 = ODEProblem(sys, u0, tend, [DOX => 295μM])
prob1 = ODEProblem(sys, u0, tend, [DOX => 295μM, ρC4 => 500μM])
prob2 = ODEProblem(sys, u0, tend, [DOX => 295μM, ρC3 => 500μM])
@time sol0 = solve(prob0, alg; opts...)
@time sol1 = solve(prob1, alg; opts...)
@time sol2 = solve(prob2, alg; opts...)

#---
plot(sol0, idxs=(s, vm), lab="DOX=295", title="PM potential") |> PNG

#---
plot(sol1, idxs=(s, vm), lab="DOX=295, ρC4=500", title="PM potential") |> PNG

#---
plot(sol2, idxs=(s, vm), lab="DOX=295, ρC3=500", title="PM potential") |> PNG

#---
fig = plot(sol0, idxs=(s, [atp_i / adp_i]), label="DOX=295", title="ATP:ADP")
plot!(fig, sol1, idxs=(s, [atp_i / adp_i]), label="C4 500uM")
plot!(fig, sol2, idxs=(s, [atp_i / adp_i]), label="C3 500uM")
fig |> PNG

#---
idxs = (sys.t/1000, sys.vO2 + sys.vROS)
fig = plot(sol0, idxs=idxs, label="DOX=295", title="Oxygen consumption")
plot!(fig, sol1, idxs=idxs, label="C4 500")
plot!(fig, sol2, idxs=idxs, label="C3 500", xlabel="Time (s)", ylabel="vO2 (μM/ms)")
fig |> PNG

# O2 shunt
idxs = (sys.t/1000, sys.vROS / (sys.vO2 + sys.vROS))
fig = plot(sol0, idxs=idxs, label="DOX=295", title="ROS shunt fraction")
plot!(fig, sol1, idxs=idxs, label="C4 500")
plot!(fig, sol2, idxs=idxs, label="C3 500", xlabel="Time (s)", legend=:right)
fig |> PNG

# Q cycle : reduced Q pool
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol0, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="Q cycle (DOX=295uM)", legend=:right) |> PNG

# Succinate accumulation and CAC slowed down.
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal= sys
plot(sol0, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], leg=:right, title="CAC metabolites") |> PNG

# NADH production from CAC decreased
idxs = (sys.t/1000, sys.nadh_m)
fig = plot(sol0, idxs=idxs, label="DOX=295", title="NADH (mito)")
plot!(fig, sol1, idxs=idxs, label="C4 500")
plot!(fig, sol2, idxs=idxs, label="C3 500", xlabel="Time (s)")
fig |> PNG
