# Faraday constant (C/mol)
const F = 96485
# Temperature (Kelvin)
const T = 310
# Ideal gas constant (kJ / (K*mol))
const R = 8314.46
# Unit voltage for Nernst potential (mV)
const RT_F = (R * T) / F
# Scaling factor of unit voltage for Nernst potential(1/mV)
const F_RT = F / (R * T)

# Dissociation constants (units in mM)
const KWATER = 1E-14 * 1E3 * 1E3
const KA_PI = 1.78E-7 * 1E3
const KA_ATP = 3.31E-7 * 1E3
const KA_ADP = 4.17E-7 * 1E3
const KMG_ATP = 6.46E-5 * 1E3
const KMG_ADP = 5.62E-4 * 1E3
const KA_SUC = 6.3E-6 * 1E3
const KA_H2O = 1E-14 * 1E3

# Hydroxide ion
_oh(h) = KWATER / h

# Michaelis-Menton reaction rates and Hills function
_mm(x, k) = x / (x + k)
_mm_reci(x, k) = (x + k) / x
_hills(x, k, n) = x^n / (x^n + k^n)
_hills_reci(x, k, n) = (x^n + k^n) / x^n

#=
Get propotions of AXP from dissociation constants
Returns AXPn-, HAXP, MgAXP, poly (the denominator)
Scalar version
=#
_poly(h, KA) = one(h) + h / KA
function _breakdown_axp(axp, h, mg, KA, KMG)
    poly = _poly(h, KA) + mg / KMG
    axpn = axp / poly
    haxp = axpn * h / KA
    mgaxp = axp - axpn - haxp
    return (axpn, haxp, mgaxp, poly)
end

breakdown_atp(atp, h, mg) = _breakdown_axp(atp, h, mg, KA_ATP, KMG_ATP)
breakdown_adp(adp, h, mg) = _breakdown_axp(adp, h, mg, KA_ADP, KMG_ADP)

# Binding polynomial of phosphate
_pi_poly(h) = _poly(h, KA_PI)
_H2PO4(pi_t, h) = pi_t * _mm(h, KA_PI)

# Binding Polynomial of succinate
_suc_poly(h) = _poly(h, KA_SUC)
# Binding Polynomial of water
_h2o_poly(h) = _poly(h, KA_H2O)

# Nernst potentials (in mV)
_nernst(x_out, x_in, z=1) = RT_F / z * log(x_out / x_in)
_ena(na_o, na_i) = _nernst(na_o, na_i, 1)
_enak(k_o, na_o, k_i, na_i, P_NA_K) = _nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i, 1)
_ek(k_o, k_i) = _nernst(k_o, k_i, 1)
_eca(ca_o, ca_i) = _nernst(ca_o, ca_i, 2)
_esox(sox_i, sox_m) = _nernst(sox_i, sox_m, -1)

# Inverse function of the Nernst equation
# ra stands for Ratio of (electrochemical) Activity
_ra(v) = exp(v * F_RT)

# Proton motive force
_ΔμH(dpsi, h_i, h_m) = dpsi + _nernst(h_i, h_m, 1)
_ΔμH(dpsi, ΔpH) = dpsi - 2.303 * RT_F * ΔpH

# Capacitance of palsma membrane (μF/cm²)
const CM = 1
const CM_INV = inv(CM)
# Capacitance of mitochondrial inner membrane (mV/ms)
const CM_MITO = 1.812E-3
const CM_MITO_INV = inv(CM_MITO)

# Simple diffusion rate
_vdiff(a, b, c_diff) = c_diff * (a - b)

# GHK flux equation
_ghk(zvfrt, ezvfrtm1, px, x_i, x_o, z=1) = px * z * zvfrt * F * ((ezvfrtm1 + 1) * x_i - x_o) / ezvfrtm1
