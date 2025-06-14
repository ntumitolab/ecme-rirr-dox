#===
# Doxorubicin addition test
===#
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DisplayAs: PNG
using ECMEDox
using ECMEDox: second, mM, Hz, μM, build_stim_callbacks
Plots.default(lw=1.5, size=(600, 600))

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl=0, tend)
@unpack DOX = sys
sts = unknowns(sys)
u0 = build_u0(sys)
alg = KenCarp47()
stim = build_stim_callbacks(sys.iStim, tend)
opts = (; reltol=1e-6, abstol=1e-6, progress=true, callback=stim)
prob = ODEProblem(sys, u0, tend)

# The collapse of MMP is between 290uM and 300uM of DOX
doxrange = 290μM:1μM:300μM
prob_func = (prob, i, repeat) -> begin
    prob.ps[DOX] = doxrange[i]
    prob
end

eprob = EnsembleProblem(prob; prob_func)
@time sim = solve(eprob, alg; trajectories=length(doxrange), opts...)

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
plot!(fig, ylabel="Conc. (μM)", xlabel="Time (s)", legend=:bottomright, ylim=(0, 8mM)) |> PNG

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
@unpack Q_n, SQn, QH2_n, QH2_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sim[end], idxs=[Q_n + Q_p, SQn, QH2_n+QH2_p], title="Q cycle ", legend=:right) |> PNG

#---
plot(sim[end], idxs=[fes_ox, fes_rd, cytc_ox, cytc_rd], title="Q cycle (downstream)", legend=:right) |> PNG
#---
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal, nadh_m = sys
plot(sim[end], idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal, nadh_m], title="TCA cycle ", legend=:right) |> PNG
