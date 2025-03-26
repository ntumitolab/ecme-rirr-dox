# # Ischemic reperfusion damage
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using DiffEqCallbacks
using ECMEDox
using ECMEDox: second, Hz, μM, nM
using Plots
using DisplayAs: PNG

#---
tend = 100.0second
bcl = 1.0second
@named sys = build_model(; bcl, tend)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
@unpack O2, J_IMAC = sys
prob = ODEProblem(sys, u0, tend, [O2 => 6nM])

reoxygen! = (integrator) -> begin
    integrator.ps[sys.O2] = 6μM
    set_proposed_dt!(integrator, 1e-6)
end

cb = PresetTimeCallback(60.0second, reoxygen!)

#---
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true, dt=1e-6, callback=cb)

#---
@unpack vm, dpsi, atp_i, adp_i = sys
pl_mmp = plot(sol, idxs=dpsi, lab=false, title="(A) Mito. memb. potential", xlabel="", ylabel="Voltage (mV)")
pl_mmp = vline!(pl_mmp, [120.0second], lines=(:dash, :black), lab=false)
pl_vm = plot(sol, idxs=vm, lab=false, title="(B) Action potential", xlabel="", ylabel="Voltage (mV)")
pl_vm = vline!(pl_vm, [120.0second], lines=(:dash, :black), lab=false)
pl_atp = plot(sol, idxs=[atp_i/ adp_i], lab=false, title="(C) ATP:ADP", xlabel="", ylabel="Ratio")
pl_atp = vline!(pl_atp, [120.0second], lines=(:dash, :black), lab=false)
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
pl_cac = plot(sol, idxs=sys.suc, legend=false, title="(D) Succinate", xlabel="Time (ms)", ylabel="Conc. (μM)")
pl_cac = vline!(pl_cac, [120.0second], lines=(:dash, :black), lab=false)
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
pl_q = plot(sol, idxs=[Q_n, Q_p, Qdot_n, QH2_n, QH2_p, Qdot_p], title="(E) Q cycle", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")
pl_q = vline!(pl_q, [120.0second], lines=(:dash, :black), lab=false)
pl_ros = plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="(F) ROS generation", lab=false, xlabel="Time (ms)", ylabel="Fraction of O2 consumption (%)")
pl_ros = vline!(pl_ros, [120.0second], lines=(:dash, :black), lab=false)
plot(pl_mmp, pl_vm, pl_atp, pl_cac, pl_q, pl_ros, size=(1200, 800)) |> PNG
