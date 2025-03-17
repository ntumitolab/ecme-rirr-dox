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
@unpack DOX, ET_SOD_I, J_IMAC = sys
sts = unknowns(sys)
u0 = build_u0(sys)
alg = KenCarp47()
opts = (; reltol=1e-6, abstol=1e-6, progress=true)

prob = ODEProblem(sys, u0, tend)

sol = solve(prob, alg; opts...)  ## Warm-up
plot(sol, idxs=[sys.sox_i, sys.sox_m])
plot(sol, idxs=[sys.vHres, sys.vHu + sys.vANT, sys.vNaCa + 2sys.vUni, sys.vIMAC])

prob300 = ODEProblem(sys, u0, tend, [DOX => 310μM, ET_SOD_I=>1μM, J_IMAC=>1.0])
sol300 = solve(prob300, alg; opts...)
sol300[sys.dpsi]

o2shunt = (sys.t/1000, sys.vROS / (sys.vO2 + sys.vROS))
plot(sol300, idxs=sys.dpsi)
plot(sol300, idxs=[sys.sox_i, sys.sox_m])
plot(sol300, idxs=o2shunt)

# IMAC cannot be closed
plot(sol300, idxs=[sys.vHres, sys.vHu, sys.vIMAC])

# SOD saturation
plot(sol300, idxs=sys.vSOD_i)

# Cell death?
plot(sol300, idxs=[sys.atp_i/sys.adp_i])
