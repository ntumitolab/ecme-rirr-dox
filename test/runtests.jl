using ECMEPCr

using DifferentialEquations

using Plots

params = CMCParams()

u0 = build_u0()

tspan = (0.0, 400.0)

cbs = build_pacing_callbacks(tspan)

prob = ODEProblem(model!, u0, tspan, params)

sol = solve(prob, Rodas5(), callback = cbs)