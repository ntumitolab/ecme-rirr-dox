using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ECMEDox
using ECMEDox: get_cicr40_sys

@parameters A_CAP V_SS_SINGLE ca_i ca_jsr ca_o vm
sys = get_cicr40_sys(ca_i, ca_jsr, ca_o, vm, A_CAP, V_SS_SINGLE)

simp = structural_simplify(sys)

observed(simp)

equations(simp)


# ASPARTATE

kfaat = 6.44E-4
oaa = 6.6423E-8
glu = 10
kasp = 1.5E-6
keq = 6.6
akg = 7.0596E-5

asp = kfaat * oaa * glu / (kasp + kfaat * akg / keq)

# SDH
g = -10.12e3
keq = exp(-g / 8.314 / 310)

asp / glu
