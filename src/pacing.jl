# The main workhorse
using Plots, DifferentialEquations
include("ecme-dox-3.jl")

u0 = collect(u0_tup)

param = Params(
  pDOX = DOXParams(DOX=0.0),
  pCK=CKParams(V_ATPASE_CYTO=1E-5),
  pC1=C1Params(SCALE=1e3 * 35),
  pC3=C3Params(SCALE=60e3 * 10),
  pC4=C4Params(SCALE=60e3 * 10),
  pC5 = C5Params(ρF1=0.5),
  pANT = ANTParams(VMAX=5E-3),
  pSODi = SODParams(ET=1.43E-3 * 0.65),
  BCL=1000.0)

@unpack A_CAP_V_MYO_F = param

tspan = (0.0, 50000.0)
prob = ODEProblem(rhs, u0, tspan, param)

integrator = init(prob, Rodas5(); reltol=1e-9, abstol=1e-9, dt=0.01, progress=true, tstops=tspan[1]:1000.0:tspan[end])
for (u,t) in tuples(integrator)
    (vm, dpsi, m_na, h_na, j_na, x_k, na_i, k_i, ca_i, ca_jsr, ca_nsr, ca_ss, ltr_ca, htr_ca,
     po1, po2, pc2, c0, c1, c2, c3, c4, o, cca0, cca1, cca2, cca3, y_ca, p0, p1, p2, p3, n1,
     adp_i, adp_ic, crp_i, crp_ic, ca_m, adp_m, nadh, isoc, akg, scoa, suc, fum, mal, oaa,
     sox_m, sox_i, h2o2_i, gsh_i, Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, b1, b2, b3, b4,
     fes_ox, cytc1_ox, cytc_ox) = u
    @show t, vm, dpsi
    @show adp_i, adp_m
    @show Q_n, QH2_n
end

solPacing = solve(prob, Rodas5(); reltol=1e-9, abstol=1e-9, dt=0.01, progress=true, tstops=tspan[1]:1000.0:tspan[end])
plot(solPacing, vars=(0, nameLUT[:vm]), label="Membrane Potential", lw=1)
plot(solPacing, vars=(0, nameLUT[:dpsi]), label="Mitochondrial Potential", lw=1)
plot(solPacing, vars=(0, [nameLUT[:adp_i], nameLUT[:adp_m]]), label=["ADPi" "ADPm"], lw=1)
plot(solPacing, vars=(0, nameLUT[:nadh]), label="[NADH]", lw=1)
plot(solPacing, vars=(0, nameLUT[:ca_m]), label="[Ca]m", lw=1)
plot(solPacing, vars=(0, nameLUT[:isoc]), label="[ISOC]", lw=1)
plot(solPacing, vars=(0, nameLUT[:akg]), label="[AKG]", lw=1)
plot(solPacing, vars=(0, nameLUT[:scoa]), label="[SCoA]", lw=1)
plot(solPacing, vars=(0, nameLUT[:suc]), label="[SUC]", lw=1)
plot(solPacing, vars=(0, nameLUT[:fum]), label="[FUM]", lw=1)
plot(solPacing, vars=(0, nameLUT[:mal]), label="[MAL]", lw=1)
plot(solPacing, vars=(0, nameLUT[:oaa]), label="[OAA]", lw=1)
plot(solPacing, vars=(0, [nameLUT[:Q_n], nameLUT[:QH2_n]]), lw=1, label=["Q_n" "QH2_n"])
plot(solPacing, vars=(0, [nameLUT[:b1], nameLUT[:b2], nameLUT[:b3], nameLUT[:b4]]), lw=1)
plot(solPacing, vars=(0, [nameLUT[:fes_ox], nameLUT[:cytc1_ox], nameLUT[:cytc_ox]]), lw=1, label=["fes_ox" "cytc1_ox" "cytc_ox"])
plot(solPacing, vars=(0, [nameLUT[:sox_i], nameLUT[:sox_m]]), lw=1, label=["sox_i" "sox_m"])
plot(solPacing, vars=(0, nameLUT[:sox_i]), lw=1, label=["sox_i"])
plot(solPacing, vars=(0, nameLUT[:gsh_i]), lw=1)
plot(solPacing, vars=(0, nameLUT[:ca_i]), lw=1)
plot(solPacing, vars=(0, nameLUT[:na_i]), lw=1)

# Mitochondrial potential going crazy!
plot(solPacing, vars=(0, nameLUT[:dpsi]), label="Mitochondrial Potential", lw=1, xlims=(0.0, 10000.0), ylims=(-10.0, 200.0))
plot(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vHu]], solPacing.t[1], solPacing.t[end], label="vHu")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vATPase]], solPacing.t[1], solPacing.t[end], label="vATPase")
plot!(solPacing, vars=(0, nameLUT[:adp_m]), lw=1)

plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:jUp]], solPacing.t[1], solPacing.t[end], label="jUp")
plot!(t-> A_CAP_V_MYO_F * ecme_dox(solPacing(t), param, t)[rateMap[:iPCa]], solPacing.t[1], solPacing.t[end], label="iPCa")
plot!(t-> A_CAP_V_MYO_F * ecme_dox(solPacing(t), param, t)[rateMap[:iNaK]], solPacing.t[1], solPacing.t[end], label="iNaK")
plot!(t-> A_CAP_V_MYO_F * ecme_dox(solPacing(t), param, t)[rateMap[:iNaCa]], solPacing.t[1], solPacing.t[end], label="iNaCa")
plot(t-> ecme_dox(solPacing(t), param, t)[rateMap[:vAM]], solPacing.t[1], solPacing.t[end], label="vAM")
plot!(t-> ecme_dox(solPacing(t), param, t)[rateMap[:d_adp_i]], solPacing.t[1], solPacing.t[end], label="d_adp_i")
