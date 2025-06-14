# # Hypoxia reperfusion simulation
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using DiffEqCallbacks
using ECMEDox
using ECMEDox: second, Hz, μM, nM, build_stim_callbacks
using Plots
using DisplayAs: PNG
Plots.default(lw=1.5)

#---
tend = 200.0second
bcl = 1.0second
@named sys = build_model(; bcl=0, tend)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
@unpack O2 = sys
prob = ODEProblem(sys, [u0;  O2 => 6nM], tend)

stim = build_stim_callbacks(sys.iStim, tend)
reoxygen! = (integrator) -> begin
    integrator.ps[sys.O2] = 6μM
    set_proposed_dt!(integrator, 1e-6)
end

cbs = CallbackSet(PresetTimeCallback(100.0second, reoxygen!), stim)

#---
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true, dt=1e-6, callback=cbs)

#---
@unpack vm, dpsi, atp_i, adp_i = sys
pl_mmp = plot(sol, idxs=dpsi, lab=false, title="(A) Mito. memb. potential", xlabel="", ylabel="Voltage (mV)")
pl_mmp = vline!(pl_mmp, [100.0second], lines=(:dash, :black), lab=false)
pl_vm = plot(sol, idxs=vm, lab=false, title="(B) Action potential", xlabel="", ylabel="Voltage (mV)")
pl_vm = vline!(pl_vm, [100.0second], lines=(:dash, :black), lab=false)
pl_atp = plot(sol, idxs=atp_i/1000, lab=false, title="(C) ATP", xlabel="", ylabel="Conc. (mM)")
pl_atp = vline!(pl_atp, [100.0second], lines=(:dash, :black), lab=false)
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
pl_cac = plot(sol, idxs=sys.suc, legend=false, title="(D) Succinate", xlabel="Time (ms)", ylabel="Conc. (μM)")
pl_cac = vline!(pl_cac, [100.0second], lines=(:dash, :black), lab=false)
@unpack Q_n, SQn, QH2_n, QH2_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
pl_q = plot(sol, idxs=[Q_n, Q_p, SQn, QH2_n, QH2_p], title="(E) Q cycle", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")
pl_q = vline!(pl_q, [100.0second], lines=(:dash, :black), lab=false)
pl_ros = plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="(F) ROS generation", lab=false, xlabel="Time (ms)", ylabel="Fraction of O2 consumption (%)")
pl_ros = vline!(pl_ros, [100.0second], lines=(:dash, :black), lab=false)
plot(pl_mmp, pl_vm, pl_atp, pl_cac, pl_q, pl_ros, size=(1200, 800)) |> PNG

#---
plot(sol, idxs=sys.nadh_m, title="NADH (mito)")

#---
plot(sol, idxs=[sys.vROSC1, sys.vROSC3])

#---
plot(sol, idxs=[sys.sox_m, sys.sox_i])
