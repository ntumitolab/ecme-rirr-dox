#===
Physical constants and common functions
===#
const ms = 1                  # millisecond
const second = 1000ms         # second is the SI unit
const minute = 60second       # minute
const metre = 1               # meter
const mmol = 1                # millimole
const μM = mmol // metre^3    # micromolar = 1
const Hz = 1 // second        # Herz
const kHz = 1 // ms           # kilohertz
const cm = metre // 100       # centimeter
const cm² = cm^2              # square centimeter
const μm = metre // 10^6      # micrometer
const mL = cm^3               # milliliter = cubic centimeter
const Liter = 1000mL          # liter
const μL = μm^3               # microliter
const pL = Liter // 10^12     # picoliter
const mol = 1000mmol          # mole
const mM = 1000μM             # mM is the SI unit
const Molar = 1000mM          # Molarity is used in equilibrium constants
const nM = μM // 1000         # nanomolar
const Ampere = 1              # current unit Ampere
const mA = Ampere // 1000     # milliampere
const μA = Ampere // 10^6     # micropampere
const Joule = 10^6            # energy unit Joule
const Kelvin = 1              # temperature unit Kelvin
const Columb = Ampere * second # unit of electric charge
const Volt = Joule // Columb   # voltage
const mV = Volt // 1000         # millivolt (mV) = 1
const milliseimens = Ampere // Volt // 1000 # milliseimens
const Farad = Columb // Volt   # Capacitance = 1
const μF = Farad // 10^6       # Microfarad
const μAμF = μA // μF           # Common unit for current density, normalized by capacitance = 1
const mSμF = milliseimens // μF # Common unit for conductance, normalized by capacitance = 1
const Faraday = 96485Columb / mol # Faraday constant (columb / mol)
const T₀ = 310Kelvin          # Default temp (37C)
const RGAS = 8.314Joule / Kelvin / mol # Ideal gas constant (J/K⋅mol)
const VT = RGAS * T₀ / Faraday # Thermal voltage (@37C), 26.7 mV
const iVT = inv(VT) # Reciprocal of thermal voltage (0.037 per motive)

# Dissociation and affinity constants
const KWATER = 1E-14 * Molar^2
const KA_PI = 1.78E-7 * Molar
const KA_ATP = 3.31E-7 * Molar
const KA_ADP = 4.17E-7 * Molar
const KMG_ATP = 6.46E-5 * Molar
const KMG_ADP = 5.62E-4 * Molar
const KA_SUC = 6.3E-6 * Molar
const KA_H2O = 1E-14 * Molar
const iKWATER = inv(KWATER)
const iKA_PI = inv(KA_PI)
const iKA_ATP = inv(KA_ATP)
const iKA_ADP = inv(KA_ADP)
const iKMG_ATP = inv(KMG_ATP)
const iKMG_ADP = inv(KMG_ADP)
const iKA_SUC = inv(KA_SUC)
const iKA_H2O = inv(KA_H2O)

#===
Commonly-used functions
===#

"""
Regular Hill/MM function
"""
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(NaNMath.pow(x, n), NaNMath.pow(k, n))

"""
Logistic sigmoid function.

    expit(x[, a=1, b=1]) = a / (b + exp(-x))
"""
expit(x, a=one(x), b=one(x)) = a / (b + exp(-x))

"""
Relative exponential function.

    exprel(x) = x / (exp(x) - 1)
"""
exprel(x) = x / expm1(x)

"Acid dissociation polynomial"
_poly(h, iKA) = h * iKA + 1
_poly(h, iKA, mg, iKMG) = h * iKA + +mg * iKMG + 1
#=
Get propotions of AXP from dissociation constants
Returns AXPn-, HAXP, MgAXP, poly (the denominator polynomial)
=#

function _breakdown_axp(axp, h, mg, iKA, iKMG)
    poly = _poly(h, iKA, mg, iKMG)
    axpn = axp / poly
    haxp = axpn * h * iKA
    mgaxp = axp - axpn - haxp
    return (axpn, haxp, mgaxp, poly)
end

breakdown_atp(atp, h, mg) = _breakdown_axp(atp, h, mg, iKA_ATP, iKMG_ATP)
breakdown_adp(adp, h, mg) = _breakdown_axp(adp, h, mg, iKA_ADP, iKMG_ADP)

"Binding polynomial of phosphate"
pipoly(h) = _poly(h, iKA_PI)

"Dihydrogen phosphate"
dhpi(Σpi, h) = Σpi * hil(h, KA_PI)

"Binding Polynomial of succinate"
sucpoly(h) = _poly(h, iKA_SUC)

"Binding Polynomial of water"
h2opoly(h) = _poly(h, iKA_H2O)

"Nernst potential"
nernst(x_out, x_in) = VT * NaNMath.log(x_out / x_in)
nernst(x_out, x_in, z) = nernst(x_out, x_in) / z
nernstNaK(k_o, na_o, k_i, na_i, P_NA_K) = nernst(na_o * P_NA_K + k_o, na_i * P_NA_K + k_i)

"Proton motive force"
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
function ghk(px, vm, x_i, x_o, z=1)
    zvfrt = z * vm * iVT
    return px * z * Faraday * exprel(zvfrt) * (exp(zvfrt) * x_i - x_o)
end
