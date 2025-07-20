# # Hypoxia reperfusion simulation
using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using DiffEqCallbacks
using ECMEDox
using ECMEDox: second, Hz, μM, nM, mV
using Plots
using DisplayAs: PNG
Plots.default(lw=1.5)

#---
tend = 100.0second
bcl = 1.0second
@named sys = build_model(; bcl, tend)
u0 = build_u0(sys)
sts = unknowns(sys)
alg = KenCarp47()
prob = ODEProblem(sys, [u0; sys.O2_o => 6nM], tend)

reoxygen! = (integrator) -> begin
    integrator.ps[sys.O2_o] = 6μM
    set_proposed_dt!(integrator, 1e-6)
end

cbs = PresetTimeCallback(50.0second, reoxygen!)

#---
@time sol = solve(prob, alg; reltol=1e-6, abstol=1e-6, progress=true, dt=1e-6, callback=cbs)

#---
@unpack vm, dpsi, atp_i, adp_i = sys
pl_mmp = plot(sol, idxs=dpsi, lab=false, title="(A) Mito. memb. potential", xlabel="", ylabel="Voltage (mV)")
pl_mmp = vline!(pl_mmp, [50.0second], lines=(:dash, :black), lab=false)
pl_vm = plot(sol, idxs=vm, lab=false, title="(B) Action potential", xlabel="", ylabel="Voltage (mV)")
pl_vm = vline!(pl_vm, [50.0second], lines=(:dash, :black), lab=false)
pl_atp = plot(sol, idxs=atp_i/1000, lab=false, title="(C) ATP", xlabel="", ylabel="Conc. (mM)")
pl_atp = vline!(pl_atp, [50.0second], lines=(:dash, :black), lab=false)
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
pl_cac = plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal], legend=:right, title="(D) CAC intermediates", xlabel="Time (ms)", ylabel="Conc. (μM)")
pl_cac = vline!(pl_cac, [50.0second], lines=(:dash, :black), lab=false)
@unpack Q_n, SQn, QH2_n, QH2_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
pl_q = plot(sol, idxs=[Q_n, Q_p, SQn, QH2_n, QH2_p], title="(E) Q cycle", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")
pl_q = vline!(pl_q, [50.0second], lines=(:dash, :black), lab=false)
pl_ros = plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="(F) ROS generation", lab=false, xlabel="Time (ms)", ylabel="Fraction of O2 consumption (%)", ylims=(0, 100))
pl_ros = vline!(pl_ros, [50.0second], lines=(:dash, :black), lab=false)
plot(pl_mmp, pl_vm, pl_atp, pl_cac, pl_q, pl_ros, size=(1200, 800)) |> PNG

# MMP
@unpack vm, dpsi, atp_i, adp_i = sys
pl_mmp = plot(sol, idxs=dpsi, lab=false, title="(A) Mito. memb. potential", xlabel="", ylabel="Voltage (mV)")
vline!(pl_mmp, [50.0second], lines=(:dash, :black), lab=false)
pl_mmp |> PNG

# Action potential
pl_vm = plot(sol, idxs=vm, lab=false, title="(B) Action potential", xlabel="", ylabel="Voltage (mV)")
vline!(pl_vm, [50.0second], lines=(:dash, :black), lab=false)
pl_vm |> PNG

# ATP
pl_atp = plot(sol, idxs=atp_i/1000, lab=false, title="(C) ATP", xlabel="", ylabel="Conc. (mM)")
vline!(pl_atp, [50.0second], lines=(:dash, :black), lab=false)
pl_atp |> PNG

# CAC
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal, nadh_m = sys
pl_cac = plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal, nadh_m], legend=:right, title="(D) CAC intermediates", xlabel="Time (ms)", ylabel="Conc. (μM)")
vline!(pl_cac, [50.0second], lines=(:dash, :black), lab=false)
pl_cac |> PNG

# CAC flux
pl = plot(sol, idxs=[sys.vSDH, sys.vAAT, sys.vMDH, sys.vIDH])
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG

# Q cycle
@unpack Q_n, SQn, QH2_n, QH2_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
pl_q = plot(sol, idxs=[Q_n, Q_p, SQn, QH2_n, QH2_p], title="(E) Q cycle", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")
vline!(pl_q, [50.0second], lines=(:dash, :black), lab=false)
pl_q |> PNG

# Downstream of Q cycle
pl = plot(sol, idxs=[fes_ox, fes_rd, cytc_ox, cytc_rd], title="Q cycle downstream", legend=:left, xlabel="Time (ms)", ylabel="Conc. (μM)")
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG

# Proton pumping
pl = plot(sol, idxs=[sys.vHresC1, sys.vHresC3, sys.vHresC4], title="Proton pumping", legend=:left, xlabel="Time (ms)", ylabel="Rate (μM/ms)")
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG

# Mito oxygen concentration
pl = plot(sol, idxs=sys.O2)
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG

# O2 shunt
pl_ros = plot(sol, idxs=100 * sys.vROS / (sys.vO2 + sys.vROS), title="(F) ROS generation", lab=false, xlabel="Time (ms)", ylabel="Fraction of O2 consumption (%)", ylims=(0, 100))
vline!(pl_ros, [50.0second], lines=(:dash, :black), lab=false)
pl_ros |> PNG

# Superoxide production
pl = plot(sol, idxs=[sys.vROSC1, sys.vROSIf, sys.vROSIq, sys.vROSC3], ylims=(0.0, 0.05))
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG

# Superoxide concentrations
pl = plot(sol, idxs=[sys.sox_m, sys.sox_i])
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG


# Semiquinone from Complex III Qo
pl = plot(sol, idxs=[sys.SQp])
vline!(pl, [50.0second], lines=(:dash, :black), lab=false)
pl |> PNG
