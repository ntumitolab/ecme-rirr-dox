#===
Physical constants and common functions
===#
@constants begin
    ms = 1                  # millisecond
    second = 1000ms         # second is the SI unit
    minute = 60second       # minute
    Hz = 1 / second         # Herz
    kHz = 1 / second        # kilohertz
    metre = 1               # meter
    cm = metre / 10^2       # centimeter
    cm² = cm^2              # square centimeter
    μm = metre / 10^6       # micrometer
    mL = cm^3               # milliliter = cubic centimeter
    Liter = 1000mL          # liter
    μL = μm^3               # microliter
    pL = Liter / 10^12      # picoliter
    mmol = 1                # millimole
    mol = 1000mmol          # mole
    μM = mmol/metre^3       # micromolar
    mM = 1000μM             # mM is the SI unit
    Molar = 1000mM          # Molarity is used in equilibrium constants
    nM = μM / 10^3          # nanomolar
    Ampere = 1              # current unit Ampere
    mA = Ampere / 10^3      # milliampere
    μA = Ampere / 10^6      # micropampere
    Joule = 10^6            # energy unit Joule
    Kelvin = 1              # temperature unit Kelvin
    Columb = Ampere * second # unit of electric charge
    Volt = Joule / Columb   # voltage
    mV = Volt / 10^3        # millivolt (mV) = 1
    milliseimens = Ampere / Volt / 10^3 # milliseimens
    Farad = Columb / Volt   # Capacitance
    μF = Farad / 10^6       # Microfarad
    T₀ = 310Kelvin          # Default temp (37C)
    Faraday = 96485Columb / mol # Faraday constant (columb / mol)
    RGAS = 8.314Joule/Kelvin/mol # Ideal gas constant (J/K⋅mol)
    VT = RGAS * T₀ / Faraday # Thermal voltage (@37C), about 26.7 mV
    iVT = Faraday / RGAS * T₀ # Reciprocal of thermal voltage
    μAμF = μA / μF           # Common unit for current density, normalized by capacitance
    mSμF = milliseimens / μF # Common unit for conductance, normalized by capacitance

    # Dissociation and affinity constants
    KWATER = 1E-14 * Molar^2
    KA_PI = 1.78E-7 * Molar
    KA_ATP = 3.31E-7 * Molar
    KA_ADP = 4.17E-7 * Molar
    KMG_ATP = 6.46E-5 * Molar
    KMG_ADP = 5.62E-4 * Molar
    KA_SUC = 6.3E-6 * Molar
    KA_H2O = 1E-14 * Molar
    iKWATER = inv(KWATER)
    iKA_PI = inv(KA_PI)
    iKA_ATP = inv(KA_ATP)
    iKA_ADP = inv(KA_ADP)
    iKMG_ATP = inv(KMG_ATP)
    iKMG_ADP = inv(KMG_ADP)
    iKA_SUC = inv(KA_SUC)
    iKA_H2O = inv(KA_H2O)
end

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

    expit(x[, a=1, b=1]) = a / (b + exp(-x))
"""
expit(x, a=one(x), b=one(x)) = a / (b + exp(-x))

"""
Relative exponential function.

    exprel(x) = x / (exp(x) - 1)
"""
exprel(x) = x / expm1(x)

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
