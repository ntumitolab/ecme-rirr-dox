using ProgressLogging
using OrdinaryDiffEq
using ModelingToolkit
using Plots
using DataFrames
using CSV
using ECMEDox
using ECMEDox: second, mM, Hz, μM
Plots.default(lw=2, size=(600, 600))

tend = 1000.0second
bcl = 1second
@named sys = build_model(; bcl, tend)
@unpack DOX, ρC4, ρC3 = sys
sts = unknowns(sys)
obs = [i.lhs for i in observed(sys)]
u0 = build_u0(sys)

# The model is very sensitive to aconitase activity
# The phase transition (of the Q cycle) is between 250uM to 260uM of DOX
prob = ODEProblem(sys, u0, tend)
alg = FBDF()
@time sol = solve(prob, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

plot(sol, idxs=sys.vm, tspan=(100second, 105second))
plot(sol, idxs=[sys.atp_m / sys.adp_m])
plot(sol, idxs=sys.dpsi)
plot(sol, idxs=sys.vC5)
plot(sol, idxs=sys.vANT)
plot(sol, idxs=sys.vHres)
plot(sol, idxs=sys.vHresC3)
plot(sol, idxs=sys.nadh_m)
plot(sol, idxs=sys.vC1)
plot(sol, idxs=[sys.vIDH, sys.vKGDH, sys.vMDH])
@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal])

plot(sol, idxs=sys.ca_m)
plot(sol, idxs=sys.ca_i)
plot(sol, idxs=[sys.vUni, sys.vNaCa])
plot(sol, idxs=sys.vUni, tspan=(1second, 5second))
plot(sol, idxs=sys.ca_i, tspan=(1second, 5second))

prob0 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC4 => 325μM])
prob1 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC4 => 500μM])
prob2 = ODEProblem(sys, u0, tend, [DOX => 260μM, ρC3 => 500μM])
alg = FBDF()
@time sol0 = solve(prob0, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)
@time sol1 = solve(prob1, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)
@time sol2 = solve(prob2, alg; reltol=1e-7, abstol=1e-7, progress=true, maxiters=1e8)

df = DataFrame(sol2)

for s in obs
    df[!, Symbol(s)] = sol2[s]
end

CSV.write("dox260-c3-500.csv", df)

@unpack atp_i, adp_i, vm, na_o, na_i, atp_i, adp_i, atp_m, adp_m, sox_i, sox_m = sys
@unpack ca_i, ca_nsr, ca_jsr, ca_ss = sys
plot(sol0, idxs=vm, label="C4 1x", title="PM potential")
plot(sol0, idxs=[atp_i / adp_i], label="C4 1x", title="ATP:ADP", density=1000)

plot(sol1, idxs=vm, label="C4 3x")
plot(sol2, idxs=vm, label="C3 3x")

p1 = plot(sol, idxs=[ca_i, ca_ss], tspan=(990, tend))
p2 = plot(sol, idxs=[ca_nsr, ca_jsr], tspan=(990, tend))
plot(p1, p2, layout=(2, 1))
@unpack IK1, IK, IKp, IKatp, INa, INaCa, ICaL = sys
plot(sol0, idxs=[IK1, IK, IKp, IKatp, ICaL], tspan=(400second, tend))

plot(sol0, idxs=[atp_i / adp_i], label="C4 1x", title="ATP:ADP")
plot!(sol1, idxs=[atp_i / adp_i], label="C4 3x", line=:dash)
plot!(sol2, idxs=[atp_i / adp_i], label="C3 3x", line=:dot)

@unpack vHres, vHu, vANT, vHleak, vNaCa, vUni, vIMAC, vROS, vTrROS, vO2 = sys
@unpack Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p, fes_ox, fes_rd, cytc_ox, cytc_rd = sys
plot(sol0, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="C4 1x")
plot(sol1, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="C4 3x")
plot(sol2, idxs=[Q_n, Qdot_n, QH2_n, QH2_p, Qdot_p, Q_p], title="C3 3x")

plot(sol0, idxs=[sys.dpsi * 1000], title="C4 1x", tspan=(700, 800))
plot!(sol1, idxs=[sys.dpsi * 1000], title="C4 3x", tspan=(700, 800))
plot!(sol2, idxs=[sys.dpsi * 1000], title="C3 3x", ylims=(100, 160), tspan=(700, 800))

plot(sol, idxs=[vHres, vHu, vANT, vHleak, vNaCa, vUni, vTrROS], ylims=(0, 5))
plot(sol, idxs=vANT, tspan=(400, 500))
plot(sol, idxs=[fes_ox / fes_rd], tspan=(400, 500))
plot(sol, idxs=[cytc_ox / cytc_rd], tspan=(400, 500))
plot(sol, idxs=[sys.sox_i, sys.h2o2_i], tspan=(740, 800))
plot(sol, idxs=[sys.vSOD_i, sys.vGPX_i, sys.vGR_i, sys.vCAT], tspan=(400, 500))
plot(sol0, idxs=[sys.vROSC1, sys.vROSC3], tspan=(700, 800), title="C4 1x")
plot(sol1, idxs=[sys.vROSC1, sys.vROSC3], tspan=(700, 800), title="C4 3x")
plot(sol2, idxs=[sys.vROSC1, sys.vROSC3], tspan=(700, 800), title="C3 3x")

plot(sol0, idxs=vO2, label="C4 1x", title="vO2")
plot!(sol1, idxs=vO2, label="C4 3x", line=:dash)
plot!(sol2, idxs=vO2, label="C3 3x", line=:dot)

plot(sol0, idxs=vROS, label="C4 1x", title="vROS")
plot!(sol1, idxs=vROS, label="C4 3x")

plot(sol, idxs=[sys.vROS / (sys.vO2 + sys.vROS)])

@unpack nad_m, nadh_m, vC1, vIDH, vKGDH, vMDH = sys
plot(sol, idxs=[nad_m / nadh_m])

@unpack cit, isoc, oaa, akg, scoa, suc, fum, mal = sys
@unpack vCS, vACO, vIDH, vKGDH, vSL, vFH, vMDH, vAAT, vSDH = sys
plot(sol, idxs=[cit, isoc, oaa, akg, scoa, suc, fum, mal])

using BenchmarkTools

x = randn(10^6)
y = rand(eltype(x))

iy = inv(y)

fd(x, y) = x ./ y
fm(x, iy) = x .* iy

@benchmark fd($x, $y)
@benchmark fm($x, $iy)
