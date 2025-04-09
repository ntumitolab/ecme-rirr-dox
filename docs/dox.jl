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
alg = KenCarp47()
opts = (; reltol=1e-6, abstol=1e-6, progress=true)

# The collapse of MMP is between 295uM and 296uM of DOX
doxrange = 294μM:1μM:297μM
prob = ODEProblem(sys, u0, tend)
prob_func = (prob, i, repeat) -> remake(prob, p=[DOX => doxrange[i]])
eprob = EnsembleProblem(prob; prob_func)

@time sim = solve(eprob, alg, EnsembleThreads(); trajectories=length(doxrange), opts...)

#---
fig = plot(title="MMP")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, sys.dpsi), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="mV", xlabel="Time (s)") |> PNG

#---
fig = plot(title="ATP")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, sys.atp_i), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="Conc. (μM)", xlabel="Time (s)", legend=:right, ylim=(0, 8mM)) |> PNG

#---
fig = plot(title="ATP synthase")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, sys.vC5), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="Rate (μM/ms)", xlabel="Time (s)", legend=:right) |> PNG

#---
fig = plot(title="Cytosolic superoxide")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, sys.sox_i), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="Conc. (μM)", xlabel="Time (s)", legend=:right) |> PNG

#---
fig = plot(title="Mitochondrial superoxide")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, sys.sox_m), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="Conc. (μM)", xlabel="Time (s)", legend=:topright) |> PNG

#---
fig = plot(title="Superoxide production")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, sys.vROS), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="Rate (μM/ms)", xlabel="Time (s)", legend=:right) |> PNG

#---
fig = plot(title="O2 shunt")
for (i, dox) in enumerate(doxrange)
    plot!(fig, sim[i], idxs=(sys.t/1000, 100*sys.vROS / (sys.vO2 + sys.vROS)), lab="DOX = $(dox) μM")
end
plot!(fig, ylabel="Percentage (%)", xlabel="Time (s)", legend=:right) |> PNG

#----
@unpack Q_n, SQn, QH2_n, QH2_p, SQp, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sim[end], idxs=[Q_n + Q_p, SQn, QH2_n+QH2_p, SQp], title="Q cycle ", legend=:right) |> PNG

#---
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
plot(sim[end], idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], title="TCA cycle ", legend=:right) |> PNG

#---
prob0 = ODEProblem(sys, u0, tend, [DOX => 296μM])
prob1 = ODEProblem(sys, u0, tend, [DOX => 296μM, ρC4 => 500μM])
prob2 = ODEProblem(sys, u0, tend, [DOX => 296μM, ρC3 => 500μM])
@time sol0 = solve(prob0, alg; opts...)
@time sol1 = solve(prob1, alg; opts...)
@time sol2 = solve(prob2, alg; opts...)

#---
s = sys.t / 1000
@unpack vm, atp_i, adp_i, dpsi = sys
plot(sol0, idxs=(s, vm), lab="DOX=296", title="PM potential") |> PNG

#---
plot(sol1, idxs=(s, vm), lab="DOX=296, ρC4=500", title="PM potential") |> PNG

#---
plot(sol2, idxs=(s, vm), lab="DOX=296, ρC3=500", title="PM potential") |> PNG

#---
i = (s, [atp_i / adp_i])
fig = plot(sol0, idxs=i, label="DOX=296", title="ATP:ADP")
plot!(fig, sol1, idxs=i, label="C4 500uM")
plot!(fig, sol2, idxs=i, label="C3 500uM")
fig |> PNG

#---
i = (s, sys.dpsi)
fig = plot(sol0, idxs=i, label="DOX=296", title="MMP")
plot!(fig, sol1, idxs=i, label="C4 500uM")
plot!(fig, sol2, idxs=i, label="C3 500uM")
fig |> PNG

#---
idxs = (sys.t/1000, sys.vO2 + sys.vROS)
fig = plot(sol0, idxs=idxs, label="DOX=296", title="Oxygen consumption")
plot!(fig, sol1, idxs=idxs, label="C4 500")
plot!(fig, sol2, idxs=idxs, label="C3 500", xlabel="Time (s)", ylabel="vO2 (μM/ms)")
fig |> PNG

# O2 shunt
idxs = (sys.t/1000, sys.vROS / (sys.vO2 + sys.vROS))
fig = plot(sol0, idxs=idxs, label="DOX=296", title="ROS shunt fraction")
plot!(fig, sol1, idxs=idxs, label="C4 500")
plot!(fig, sol2, idxs=idxs, label="C3 500", xlabel="Time (s)", legend=:right)
fig |> PNG

# Q cycle : reduced Q pool (QH2 accumulation)
@unpack Q_n, SQn, QH2_n, QH2_p, SQp, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol0, idxs=[Q_n + Q_p, SQn, QH2_n+QH2_p, SQp], title="Q cycle (DOX=296uM)", legend=:right) |> PNG

# Succinate accumulation and CAC slowed down.
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal= sys
plot(sol0, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], legend=:right, title="CAC metabolites") |> PNG

# TCA cycle flux
@unpack vCS, vACO, vIDH, vKGDH, vSL, vFH, vMDH, vAAT, vSDH = sys
plot(sol0, idxs=[vCS, vACO, vIDH, vKGDH, vFH, vMDH, vAAT, vSDH], legend=:right, title="CAC flux", tspan=(100second, 1000second)) |> PNG

# ROS production from complex I and complex III
plot(sol0, idxs=[sys.vROSC1, sys.vROSC3], legend=:right, title="Superoxide generation") |> PNG

# Reverse electron transport in complex I?
plot(sol0, idxs=[sys.vQC1], legend=:right, title="Complex I rates") |> PNG

# NADH production from CAC decreased
idxs = (sys.t/1000, sys.nadh_m)
fig = plot(sol0, idxs=idxs, label="DOX=296", title="NADH (mito)")
plot!(fig, sol1, idxs=idxs, label="C4 500")
plot!(fig, sol2, idxs=idxs, label="C3 500", xlabel="Time (s)", legend=:right)
fig |> PNG
