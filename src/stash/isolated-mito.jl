# The main workhorse
using Plots, DifferentialEquations, StaticArrays
using JLD2
pyplot()

include("ecme-dox.jl")
u0 = collect(u0_tup)


function isolated_model(u, p, t)
  (d_vm, d_dpsi, d_m_na, d_h_na, d_j_na, d_x_k, d_na_i, d_k_i, d_ca_i, d_ca_jsr, d_ca_nsr, d_ca_ss, d_ltr_ca, d_htr_ca,
     d_po1, d_po2, d_pc2, d_c0, d_c1, d_c2, d_c3, d_c4, d_o, d_cca0, d_cca1, d_cca2, d_cca3, d_y_ca, d_p0, d_p1, d_p2, d_p3, d_n1,
     d_adp_i, d_adp_ic, d_crp_i, d_crp_ic, d_ca_m, d_adp_m, d_nadh, d_isoc, d_akg, d_scoa, d_suc, d_fum, d_mal, d_oaa,
     d_sox_m, d_sox_i, d_h2o2_i, d_gsh_i, d_Q_n, d_Qdot_n, d_QH2_n, d_QH2_p, d_Qdot_p, d_Q_p, d_b1, d_b2, d_b3, d_b4,
     d_fes_ox, d_cytc1_ox, d_cytc_ox) = rhs(u, p, t)

     # Set cytosolic components 0
     d_vm = 0
   	 d_m_na = d_h_na = d_j_na = d_na_i = 0
   	 d_x_k = d_k_i = 0
   	 d_p0 = d_p1 = d_p2 = d_p3 = d_n1 = 0
   	 d_ca_i = d_ca_jsr = d_ca_nsr = d_ca_ss = d_ltr_ca = d_htr_ca = 0
   	 d_adp_i = d_adp_ic = d_crp_i = d_crp_ic = 0
     d_po1 = d_po2 = d_pc2 = d_c0 = d_c1 = d_c2 = d_c3 = d_c4 = d_o = d_cca0 = d_cca1 = d_cca2 = d_cca3 = d_y_ca = 0
     return @SVector[d_vm, d_dpsi, d_m_na, d_h_na, d_j_na, d_x_k, d_na_i, d_k_i, d_ca_i, d_ca_jsr, d_ca_nsr, d_ca_ss, d_ltr_ca, d_htr_ca,
        d_po1, d_po2, d_pc2, d_c0, d_c1, d_c2, d_c3, d_c4, d_o, d_cca0, d_cca1, d_cca2, d_cca3, d_y_ca, d_p0, d_p1, d_p2, d_p3, d_n1,
        d_adp_i, d_adp_ic, d_crp_i, d_crp_ic, d_ca_m, d_adp_m, d_nadh, d_isoc, d_akg, d_scoa, d_suc, d_fum, d_mal, d_oaa,
        d_sox_m, d_sox_i, d_h2o2_i, d_gsh_i, d_Q_n, d_Qdot_n, d_QH2_n, d_QH2_p, d_Qdot_p, d_Q_p, d_b1, d_b2, d_b3, d_b4,
        d_fes_ox, d_cytc1_ox, d_cytc_ox]
end

param = Params(
  pDOX = DOXParams(DOX=0.0),
  pCK = CKParams(V_ATPASE_CYTO= 1E-5),
  pC1 = C1Params(SCALE=1e3 * 60),
  pC3 = C3Params(SCALE=60e3 * 40),
  pC4 = C4Params(SCALE=60e3 * 30),
  pC5 = C5Params(ρF1=0.05),
  pANT = ANTParams(VMAX=5E-3),
  pSODi = SODParams(ET=1.43E-3 * 0.65),
  BCL=0)

# Run 1500 seconds
tspan = (0.0, 1.5E6)
prob = ODEProblem(isolated_model, u0, tspan, param)
baselineSol = solve(prob, Rodas5(); reltol=1e-9, abstol=1e-9, dt=0.01, progress=true, dtmax=3000.0)

paramDOX = Params(
  pDOX = DOXParams(DOX=0.25),
  pCK = CKParams(V_ATPASE_CYTO= 1E-5),
  pC1 = C1Params(SCALE=1e3 * 60),
  pC3 = C3Params(SCALE=60e3 * 40),
  pC4 = C4Params(SCALE=60e3 * 30),
  pC5 = C5Params(ρF1=0.05),
  pANT = ANTParams(VMAX=5E-3),
  pSODi = SODParams(ET=1.43E-3 * 0.65),
  BCL=0)

prob = ODEProblem(isolated_model, u0, tspan, paramDOX)
doxSol = solve(prob, Rodas5(); reltol=1e-9, abstol=1e-9, dt=0.01, progress=true, dtmax=3000.0)

JLD2.@load "mito.jld2" doxSol

pDpsi = plot(doxSol, vars=(0, nameLUT[:dpsi]), label="Mitochondrial Potential", lw=1, legend=:false, xaxis = ("time (ms)"), yaxis=("Voltage (mV)"))
savefig(pDpsi, "dpsi.png")

pATPase = plot(t-> ecme_dox(doxSol(t), paramDOX, t)[rateMap[:vATPase]], doxSol.t[1], doxSol.t[end],
     label="vATPase", ylims=(-0.003, 0.0005), xaxis = ("time (ms)"), yaxis=("ATP synthase rate (mM/ms)"), legend=:false)
savefig(pATPase, "atpase.png")

pROS = plot(doxSol, vars=(0, [nameLUT[:sox_i], nameLUT[:sox_m]]), lw=1, label=["superoxide_cyto" "superoxide_mito"], xaxis = ("time (ms)"), yaxis=("Concentration (mM)"))
savefig(pROS, "ros-burst.png")
plot!(pROS, doxSol, vars=(0, [nameLUT[:Q_n], nameLUT[:QH2_n], nameLUT[:Qdot_p], nameLUT[:Q_p], nameLUT[:QH2_p]]), lw=1, label=["Q_n" "QH2_n" "Qdot_p" "Q_p" "QH2_p"], legend=:right, xaxis = ("time (ms)"))
savefig(pROS, "ros-burst-with-Q.png")

plot()
plot(t-> ecme_dox(doxSol(t), param, t)[rateMap[:vc1]], doxSol.t[1], doxSol.t[end], label="vc1")
plot!(t-> ecme_dox(doxSol(t), param, t)[rateMap[:vROSC1]], doxSol.t[1], doxSol.t[end], label="vROSC1")
plot!(t-> ecme_dox(doxSol(t), param, t)[rateMap[:vROSC3]], doxSol.t[1], doxSol.t[end], label="vROSC3")
ylims!(-1.1e-5, 3e-5)
savefig("vROS.png")

using JLD2
@save "mito.jld2" doxSol

# Other variables you could see
plot(doxSol, vars=(0, [nameLUT[:Q_n], nameLUT[:QH2_n]]), lw=1, label=["Q_n" "QH2_n"])
plot(doxSol, vars=(0, nameLUT[:nadh]), label="[NADH]", lw=1, legend=:bottomleft)
plot(doxSol, vars=(0, nameLUT[:suc]), label="[SUC]", lw=1, legend=:topleft)
plot!(t-> ecme_dox(doxSol(t), paramDOX, t)[rateMap[:vANT]], doxSol.t[1], doxSol.t[end], label="vANT")
plot(doxSol, vars=(0, nameLUT[:adp_m]), label="[ADP]m", lw=1)
plot(doxSol, vars=(0, [nameLUT[:Q_n], nameLUT[:QH2_n]]), lw=1, label=["Q_n" "QH2_n"])
plot(doxSol, vars=(0, [nameLUT[:b1], nameLUT[:b2], nameLUT[:b3], nameLUT[:b4]]), lw=1)
plot(doxSol, vars=(0, [nameLUT[:fes_ox], nameLUT[:cytc1_ox], nameLUT[:cytc_ox]]), lw=1, label=["fes_ox" "cytc1_ox" "cytc_ox"])
plot(doxSol, vars=(0, [nameLUT[:sox_i], nameLUT[:sox_m]]), lw=1, label=["sox_i" "sox_m"])
plot(doxSol, vars=(0, nameLUT[:sox_i]), lw=1, label=["sox_i"])
plot(doxSol, vars=(0, nameLUT[:gsh_i]), lw=1)
