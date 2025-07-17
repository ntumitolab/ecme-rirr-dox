using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
using DiffEqCallbacks

ph = 7:0.1:9

h_m = exp10.(-ph)

KH1_MDH=1.131E-8
KH2_MDH=26.7E-2
KH3_MDH=6.68E-12
KH4_MDH=5.62E-9
K_OFFSET_MDH = 0.04

f_ha = @. K_OFFSET_MDH + ( KH1_MDH * KH2_MDH/ (KH1_MDH * KH2_MDH + KH2_MDH * h_m + h_m * h_m))

f_hi = @. 1 + KH3_MDH / h_m * (1 + KH4_MDH / h_m)
f = f_ha .* f_hi

plot(ph, [f_ha f_hi f], labels=["f_ha" "f_hi" "f"])
