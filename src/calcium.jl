# Good old CICR

using Parameters
import .Utils: μA, cm², mM, kHz, mV, cm, Hz, one_m, expit, ghkVm, mmr

"Calcium buffering parameters"
@with_kw struct CaBuffer
    KM_CA               # Half-saturate Calcium concentration (mM)
    ET                  # Ca buffer amount (mM)
end

# Calcium buffering factor
β_ca(ca, KM, ET) = mm((ca + KM)^2, KM * ET)
β_ca(ca, p::CaBuffer) = β_ca(ca, p.KM_CA, p.ET)
(p::CaBuffer)(ca) = β_ca(ca, p)

"L-type calcium channel (LCC) parameters"
@with_kw struct LCC
    A = 2.0
    B = 2.0
    ω = 0.01kHz
    f = 0.3kHz
    g = 2.0kHz
    P_CA = 1.24E-3cm * Hz
    P_K = 1.11E-11cm * Hz
    I_CA_HALF = -0.4583μA / cm²
end

"LCC ODE system"
function lcc_system(c0, c1, c2, c3, c4, o, cca0, cca1, cca2, cca3, ca_ss, vm, p::LCC)
    @unpack A, B, ω, f, g = p
    cca4 = one_m(c0 + c1 + c2 + c3 + c4 + o + cca0 + cca1 + cca2 + cca3)
    v = vm / mV
    α = 0.4kHz * exp((v + 2) / 10)
    β = 0.05kHz * exp(-(v + 2) / 13)
    α′ = α * A
    β′ = β / B
    γ = 0.1875 * ca_ss
    v01 = 4 * α * c0 - β * c1
    v12 = 3 * α * c1 - 2 * β * c2
    v23 = 2 * α * c2 - 3 * β * c3
    v34 = α * c3 - 4 * β * c4
    v45 = f * c4 - g * o
    v67 = 4 * α′ * cca0 - β′ * cca1
    v78 = 3 * α′ * cca1 - 2 * β′ * cca2
    v89 = 2 * α′ * cca2 - 3 * β′ * cca3
    v910 = α′ * cca3 - 4 * β′ * cca4
    v06 = γ * c0 - ω * cca0
    v17 = A * γ * c1 - ω / B * cca1
    v28 = (A^2) * γ * c2 - ω / (B^2) * cca2
    v39 = (A^3) * γ * c3 - ω / (B^3) * cca3
    v410 = (A^4) * γ * c4 - ω / (B^4) * cca4
    d_c0 = -v01 - v06
    d_c1 = v01 - v12 - v17
    d_c2 = v12 - v23 - v28
    d_c3 = v23 - v34 - v39
    d_c4 = v34 - v45 - v410
    d_o = v45
    d_cca0 = v06 - v67
    d_cca1 = v17 + v67 - v78
    d_cca2 = v28 + v78 - v89
    d_cca3 = v39 + v89 - v910
    # d_cca4 = v410 + v910
    return (; d_c0, d_c1, d_c2, d_c3, d_c4, d_o, d_cca0, d_cca1, d_cca2, d_cca3)
end

"Ryanodine receptor"
@with_kw struct RyR
    R_RYR = 3.6kHz
    N = 4
    M = 3
    KA_P = 1.125E10kHz / mM^4
    KA_M = 0.576kHz
    KB_P = 4.05E6kHz / mM^3
    KB_M = 1.93kHz
    KC_P = 0.1kHz
    KC_M = 8E-4kHz
end

# Voltage inactivated factor of LCC
function d_yca(y_ca, vm)
    v = vm / mV
    y_inf = expit(-(v + 55) / 7.5) + 0.5 * expit((v - 21) / 6)
    τ = 20.0ms + 600.0ms * expit(-(v + 30) / 9.5)
    return (y_inf - y_ca) / τ
end

function i_cal(y_ca, o, vm, k_i, k_o, ca_o, p::LCC)
    @unpack P_CA, P_K, I_CA_HALF = p
    iCaMax = ghkVm(P_CA, vm, 0.001, 0.341 * ca_o, 2)
    iCaL = 6 * y_ca * o * iCaMax
    iCaK = mmr(iCaMax, I_CA_HALF) * y_ca * o * ghkVm(P_K, vm, k_i, k_o, 1)
    return (; iCaL, iCaK)
end

function ryr_system(po1, po2, pc2, ca_ss, p::RyR, CA_SS_THR = 0.095)
    @unpack N, M, KA_P, KA_M, KB_P, KB_M, KC_P, KC_M = p
    pc1 = one_m(po1 + po2 + pc2)
    # Switch formulation (rapid equlibrium assumption) by Plank et al. (2008)
    vc1o1 = -KA_M * po1 + KA_P * ca_ss^N * pc1
    vo1o2 = KB_P * ca_ss^M * po1 - KB_M * po2
    vo1c2 = KC_P * po1 - KC_M * pc2
    d_po1 = vc1o1 - vo1o2 - vo1c2
    d_po2 = vo1o2
    d_pc2 = vo1c2
    return (; d_po1, d_po2, d_pc2)
end

j_rel(po1, po2, ca_jsr, ca_ss, R_RYR) = R_RYR * (po1 + po2) * (ca_jsr - ca_ss)
j_rel(po1, po2, ca_jsr, ca_ss, p::RyR) = j_rel(po1, po2, ca_jsr, ca_ss, p.R_RYR)