#===
Physical constants and common functions
===#

# Units
const second = 1           # second
const minute = 60second    # minute
const ms = 1e-3second      # millisecond
const Hz = inv(second)     # Herz
const kHz = 1e3Hz          # kilohertz
const metre = 1            # meter
const cm = 0.01metre       # centimeter
const cm² = cm^2           # square centimeter
const μm = 1E-6metre       # Micrometer
const mL = cm^3            # milliliter = cubic centimeter
const Liter = 1e3mL        # liter
const μL = 1E-6Liter
const pL = 1E-12Liter
const mM = 1
const Molar = 1000mM       # molar (1000 since the SI units is mM)
const μM = 1E-3mM          # micromolar
const nM = 1E-6mM          # nanomolar
const Amp = 1              # ampere
const mA = 1E-3Amp         # milliampere
const μA = 1E-6Amp         # micrpampere
const Volt = 1             # volt
const mV = 1E-3Volt        # millivolt
const mS = mA / Volt       # milliseimens
const T₀ = 310.0           # Default temp (37C)
const Faraday = 96485      # Faraday constant (columb / mol)
const Farad = Amp * second / Volt
const μF = 1E-6Farad
const R = 8.314            # Ideal gas constant
const VT = R * T₀ / Faraday# Thermal voltage (@37C), around 26.7 mV
const iVT = inv(VT)        # Reciprocal of thermal voltage (@37C)

# Dissociation constants
const KWATER = 1E-14 * Molar^2
const KA_PI = 1.78E-7 * Molar
const KA_ATP = 3.31E-7 * Molar
const KA_ADP = 4.17E-7 * Molar
const KMG_ATP = 6.46E-5 * Molar
const KMG_ADP = 5.62E-4 * Molar
const KA_SUC = 6.3E-6 * Molar
const KA_H2O = 1E-14 * Molar

# Affinity constants
iKWATER = inv(KWATER)
iKA_PI = inv(KA_PI)
iKA_ATP = inv(KA_ATP)
iKA_ADP = inv(KA_ADP)
iKMG_ATP = inv(KMG_ATP)
iKMG_ADP = inv(KMG_ADP)
iKA_SUC = inv(KA_SUC)
iKA_H2O = inv(KA_H2O)

#===
Commonly-used functions
===#

"""
Regular Hill/MM function
"""
hil(x, k = one(x)) = x / (x + k)
hil(x, k, n) = hil(NaNMath.pow(x, n), NaNMath.pow(k, n))

"""
Repressive Hill/MM function
"""
hilr(x, k = one(x)) = hil(k, x)
hilr(x, k, n) = hil(k, x, n)

"""
Logistic sigmoid function.
"""
expit(x) = LogExpFunctions.logistic(x)

"Relative exponential function"
exprel(x) = x / expm1(x)
exprel2(x) = log1pexp(-log(2) * x) * inv(log(2))

#=
Get propotions of AXP from dissociation constants
Returns AXPn-, HAXP, MgAXP, poly (the denominator polynomial)
=#
_poly(h⁺, iKA) = h⁺ * iKA + 1
function _breakdown_axp(axp, h, mg, iKA, iKMG)
    poly = h⁺ * iKA + mg * iKMG + 1
    axpn = axp / poly
    haxp = axpn * h * iKA
    mgaxp = axp - axpn - haxp
    return (axpn, haxp, mgaxp, poly)
end

breakdown_atp(atp, h, mg) = _breakdown_axp(atp, h, mg, iKA_ATP, iKMG_ATP)
breakdown_adp(adp, h, mg) = _breakdown_axp(adp, h, mg, iKA_ADP, iKMG_ADP)

# Binding polynomial of phosphate
pipoly(h) = 1 + h * iKA_PI
H₂PO₄(Σpi, h) = Σpi * hil(h, KA_PI)

# Binding Polynomial of succinate
sucpoly(h) = 1 + h * iKA_SUC
# Binding Polynomial of water
h2opoly(h) = 1 + h * iKA_H2O

# Nernst potentials
nernst(x_out, x_in) = VT * NaNMath.log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) / z
nernstNaK(k_o, na_o, k_i, na_i, P_NA_K) = nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i)
Δp(ΔΨ, h_m, h_i) = ΔΨ + nernst(h_i, h_m)
Δp(ΔΨ, ΔpH) = ΔΨ - VT * log(10) * ΔpH

"""
GHK flux equation

    ghk(px, vm, x_i, x_o, z = 1)

https://en.wikipedia.org/wiki/Goldman%E2%80%93Hodgkin%E2%80%93Katz_flux_equation

- px: the permeability of the membrane for ion x measured in m·Hz
- vm: the transmembrane potential in volts
- x_i: the intracellular concentration of ion (mM)
- x_o: the extracellular concentration of ion (mM)
- z: the valence of ion x
"""
function ghk(px, vm, x_i, x_o, z = 1)
    zvfrt = z * vm * iVT
    return px * z * Faraday * exprel(zvfrt) * (exp(zvfrt) * x_i - x_o)
end

"GHK flux equation"
ghkVm(px, vm, x_i, x_o, z = 1) = ghk(px, vm, x_i, x_o, z)
