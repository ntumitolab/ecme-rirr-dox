using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ECMEDox
using ECMEDox: get_cicr40_sys

@parameters A_CAP V_SS_SINGLE ca_i ca_jsr ca_o vm
sys = get_cicr40_sys(ca_i, ca_jsr, ca_o, vm, A_CAP, V_SS_SINGLE)

simp = structural_simplify(sys)

observed(simp)

equations(simp)
