module Utils

# Units
const second = float(1)    # second
const minute = 60second    # minute
const ms = 1e-3second      # millisecond
const Hz = inv(second)     # Herz
const kHz = 1e3Hz          # kilohertz
const metre = float(1)     # meter
const cm = 0.01metre       # centimeter
const cm² = cm^2           # square centimeter
const μm = 1E-6metre       # Micrometer
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1e3mL        # liter
const μL = 1E-6Liter
const pL = 1E-12Liter
const mM = float(1)
const Molar = 1000mM       # molar (1000 since the SI units is mM)
const μM = 1E-3mM          # micromolar
const nM = 1E-6mM          # nanomolar
const Amp = float(1)       # ampere
const mA = 1E-3Amp         # milliampere
const μA = 1E-6Amp         # micrpampere
const Volt = float(1)      # volt
const mV = 1E-3Volt        # millivolt
const mS = mA / Volt       # milliseimens
const T₀ = 310.0           # Default temp (37C)
const F = 96485.0          # Faraday constant (columb / mol)
const R = 8.314            # Ideal gas constant
const VT = R * T₀ / F      # Thermal voltage (@37C), around 26.7 mV
const iVT = inv(VT)        # Reciprocal of thermal voltage (@37C)

# Capacitance of plasma membrane (μF/cm²)
const CM = 1
const iCM = inv(CM)
# Capacitance of mitochondrial inner membrane
const CM_MITO = 1.812E-3mM / mV
const iCM_MITO = inv(CM_MITO)

# Cellular Parameters
const A_CAP = 1.534E-4cm²
const V_MYO = 25.84pL
const V_MT = 15.89pL
const V_NSR = 1.4pL
const V_JSR = 0.16pL
const V_SS = 0.495E-3pL
const R_VMT_VMYO = V_MT / V_MYO           # Ratio of mitochondrial matrix to cytosolic ion diffusion space
const A_CAP_V_MYO_F = A_CAP / (V_MYO * F) # conversion factor from current (μA/cm²) to traslocation rate (mM/ms)
const A_CAP_V_SS_F = A_CAP / (V_SS * F)

# Dissociation constants
const KWATER = 1E-14 * Molar^2
const KA_PI = 1.78E-7 * Molar
const KA_ATP = 3.31E-7 * Molar
const KA_ADP = 4.17E-7 * Molar
const KMG_ATP = 6.46E-5 * Molar
const KMG_ADP = 5.62E-4 * Molar
const KA_SUC = 6.3E-6 * Molar
const KA_H2O = 1E-14 * Molar

# Incerse of Dissociation constants (affinity)
const iKWATER = inv(KWATER)
const iKA_PI = inv(KA_PI)
const iKA_ATP = inv(KA_ATP)
const iKA_ADP = inv(KA_ADP)
const iKMG_ATP = inv(KMG_ATP)
const iKMG_ADP = inv(KMG_ADP)
const iKA_SUC = inv(KA_SUC)
const iKA_H2O = inv(KA_H2O)

##################################
### Commonly-used functions
##################################

"Return `1+x` with the same type as x"
p_one(x) = one(x) + x

"Return `1-x` with the same type as x"
one_m(x) = one(x) - x

"""
Regular Hill function
"""
hill(x, k = one(x)) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)

"""
Repressive Hill function
"""
hillr(x, k = one(x)) = hill(k, x)
hillr(x, k, n) = hill(k, x, n)

"""
Logistic sigmoid function.
See scipy example https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.expit.html
"""
expit(x) = hillr(exp(-x))

"""
    exprel(x, em1 = expm1(x))

Returns `x / (exp(x)-1)` accurately when x is near zero.
See scipy example https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exprel.html
Note the fraction is the inverse of `scipy.exprel()`
"""
function exprel(x, em1 = expm1(x))
    res = x / em1
    return ifelse(x ≈ zero(x), one(res), res)
end

"""
Signed sqare root
"""
sqrt_s(x) = flipsign(sqrt(abs(x)), x)

"""
Signed power
"""
pow_s(x, n) = flipsign(abs(x)^n, x)

"""
Signed Hill function
"""
hill_s(x, k, n) = hill(pow_s(x, n), pow_s(k, n))


#=
Get propotions of AXP from dissociation constants
Returns AXPn-, HAXP, MgAXP, poly (the denominator polynomial)
=#
_poly(h⁺, iKA) = p_one(h⁺ * iKA)
function _breakdown_axp(axp, h, mg, iKA, iKMG)
    poly = p_one(h⁺ * iKA + mg * iKMG)
    axpn = axp / poly
    haxp = axpn * h * iKA
    mgaxp = axp - axpn - haxp
    return (axpn, haxp, mgaxp, poly)
end

breakdown_atp(atp, h, mg) = _breakdown_axp(atp, h, mg, iKA_ATP, iKMG_ATP)
breakdown_adp(adp, h, mg) = _breakdown_axp(adp, h, mg, iKA_ADP, iKMG_ADP)

# Binding polynomial of phosphate
pipoly(h⁺) = p_one(h⁺ * iKA_PI)
H₂PO₄(Σpi, h⁺) = Σpi * mm(h⁺, KA_PI)

# Binding Polynomial of succinate
sucpoly(h⁺) = p_one(h⁺ * iKA_SUC)
# Binding Polynomial of water
h2opoly(h⁺) = p_one(h⁺ * iKA_H2O)

# Nernst potentials
nernst(x_out, x_in) = VT * log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) / z
nernstNaK(k_o, na_o, k_i, na_i, P_NA_K) = nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i)
Δp(ΔΨ, h_m, h_i) = ΔΨ + nernst(h_i, h_m)
Δp(ΔΨ, ΔpH) = ΔΨ - VT * log(10) * ΔpH

# GHK flux equation
ghk(px, x_i, x_o, zvfrt, ezvfrtm1 = expm1(zvfrt), z = 1) = px * z * F * (p_one(ezvfrtm1) * x_i - x_o) * exprel(zvfrt, ezvfrtm1)

# GHK flux equation from voltage across the membrane
function ghkVm(px, vm, x_i, x_o, z = 1)
    zvfrt = z * vm * iVT
    em1 = expm1(zvfrt)
    return ghk(px, x_i, x_o, zvfrt, em1, z)
end

end # module Utils