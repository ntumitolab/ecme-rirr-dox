# Good old CICR with

using Parameters
include("common.jl")

@with_kw struct LCCParams
    A = 2.0
    B = 2.0
    ω = 0.01
    f = 0.3
    g = 2.0
    P_CA = 1.24E-3
    P_K = 1.11E-11
    I_CA_HALF = -0.4583
end

@with_kw struct RYRParams
    R_RYR = 3.6
    N = 4
    M = 3
    KA_P = 1.125E10
    KA_M = 0.576
    KB_P = 4.05E6
    KB_M = 1.93
    KC_P = 0.1
    KC_M = 8E-4
end

# LCC system
function lcc_system(c0, c1, c2, c3, c4, o, cca0, cca1, cca2, cca3, ca_ss, vm, pLCC::LCCParams)
    @unpack A, B, ω, f, g = pLCC
    cca4 = one(c0) - c0 - c1 - c2 - c3 - c4 - o - cca0 - cca1 - cca2 - cca3
    α = 0.4 * exp((vm + 2) / 10)
    β = 0.05 * exp(-(vm+2) / 13)
    α′ = A * α
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
    v28 = (A^2) * γ * c2 - ω / (B^2)  * cca2
    v39 = (A^3) * γ * c3 - ω  / (B^3) * cca3
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
    return (d_c0, d_c1, d_c2, d_c3, d_c4, d_o, d_cca0, d_cca1, d_cca2, d_cca3)
end

# Voltage inactivated factor of LCC
function dyca(y_ca, vm)
    y_inf = inv(1 + exp((vm+55) / 7.5)) + 0.5 / (1 + exp((-vm + 21) / 6))
    τ = 20.0 + 600.0 / (1 + exp((vm + 30) / 9.5))
    return (y_inf - y_ca) / τ
end

function ical(y_ca, o, vfrt, evfrtm1, k_i, k_o, ca_o, pLCC::LCCParams)
    @unpack P_CA, P_K, I_CA_HALF = pLCC
    ev2frtm1 = evfrtm1 * (evfrtm1 + 2)
    iCaMax = _ghk(2 * vfrt, ev2frtm1, P_CA, 0.001, 0.341 * ca_o, 2)
    iCaL = 6 * y_ca * o * iCaMax
    iCaK = _mm(I_CA_HALF, iCaMax) * y_ca * o * _ghk(vfrt, evfrtm1, P_K, k_i, k_o)
    return (iCaL, iCaK)
end

function ryr_system(po1, po2, pc2, ca_ss, pRYR::RYRParams, CA_SS_THR=0.095)
    @unpack N, M, KA_P, KA_M, KB_P, KB_M, KC_P, KC_M = pRYR
    pc1 = one(po1) - (po1 + po2 + pc2)
    # Switch formulation (rapid equlibrium assumption) by Plank et al. (2008)
    KA_P_CA = KA_P * ca_ss^N
    if ca_ss >= CA_SS_THR
        po1 = (po1 + pc1) * _mm(KA_P_CA, KA_M)
        d_po1 = zero(po1)
    else
        d_po1 = -KA_M * po1 + KA_P_CA * pc1
    end
    vo1o2 = KB_P * ca_ss^M * po1 - KB_M * po2
    vo1c2 = KC_P * po1 - KC_M * pc2
    d_po1 += -vo1o2 - vo1c2
    d_po2 = vo1o2
    d_pc2 = vo1c2
    return (d_po1, d_po2, d_pc2)
end

jrel(po1, po2, ca_jsr, ca_ss, R_RYR) = R_RYR * (po1 + po2) * (ca_jsr - ca_ss)
