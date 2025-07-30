using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ECMEDox: mV, nM, mM, kHz, ms, cm, Hz, μAcm⁻², μM, add_rate!, expit, ghk, hil
using Plots

function get_cicr16_eqs(; vm=-80mV, ca_ss=100nM, ca_jsr=1mM, k_i=144mM, ca_o=2mM, k_o=5.4mM)
    @independent_variables t
    D = Differential(t)
    @parameters begin
        A_LCC = 2
        B_LCC = 1/2
        ω_LCC = 0.01kHz
        f_LCC = 0.3kHz
        g_LCC = 2kHz
        γ_LCC = 0.1875 / (ms * mM)
        P_CA_LCC = 1E-3cm * Hz # 1.24E-3cm * Hz
        P_K_LCC = 1.11E-11cm * Hz
        I_CA_HALF_LCC = -0.4583μAcm⁻²
        R_RYR = 3.6kHz
        KA_P_RYR = 1.125E10 / (mM^4 * ms)
        KA_M_RYR = 0.576kHz
        KB_P_RYR = 4.05E6 / (mM^3 * ms)
        KB_M_RYR = 1.93kHz
        KC_P_RYR = 0.1kHz
        KC_M_RYR = 8E-4kHz
    end

    @variables begin
        # State transition rates
        α_lcc(t)
        β_lcc(t)
        # Closed LCC state
        c0_lcc(t) # 0.9991
        c1_lcc(t) = 0
        c2_lcc(t) = 0
        c3_lcc(t) = 0
        c4_lcc(t) = 0
        # Opened LCC state
        o_lcc(t) = 0
        # Ca-inhibited LCC states
        cca0_lcc(t) = 0
        cca1_lcc(t) = 0
        cca2_lcc(t) = 0
        cca3_lcc(t) = 0
        cca4_lcc(t) = 0
        # Voltage-gated LCC
        x_yca(t) = 1
        y_inf(t)
        τ_yca(t)
        # RyR
        po1_ryr(t) = 0
        po2_ryr(t) = 0
        pc1_ryr(t) ## Conserved
        pc2_ryr(t) = 0
        # Currents
        ICaMax(t)
        ICaK(t)
        ICaL(t)
        Jrel(t)
    end

    v = vm / mV
    ω, f, g = ω_LCC, f_LCC, g_LCC
    α = α_lcc
    β = β_lcc
    α′ = A_LCC * α
    β′ = B_LCC * β
    γ = γ_LCC * ca_ss
    rs_lcc = Dict()
    add_rate!(rs_lcc, 4α, c0_lcc, β, c1_lcc)
    add_rate!(rs_lcc, 3α, c1_lcc, 2β, c2_lcc)
    add_rate!(rs_lcc, 2α, c2_lcc, 3β, c3_lcc)
    add_rate!(rs_lcc, α, c3_lcc, 4β, c4_lcc)
    add_rate!(rs_lcc, f, c4_lcc, g, o_lcc)
    add_rate!(rs_lcc, 4α′, cca0_lcc, β′, cca1_lcc)
    add_rate!(rs_lcc, 3α′, cca1_lcc, 2β′, cca2_lcc)
    add_rate!(rs_lcc, 2α′, cca2_lcc, 3β′, cca3_lcc)
    add_rate!(rs_lcc, α′, cca3_lcc, 4β′, cca4_lcc)
    add_rate!(rs_lcc, γ, c0_lcc, ω, cca0_lcc)
    add_rate!(rs_lcc, A_LCC * γ, c1_lcc, B_LCC * ω, cca1_lcc)
    add_rate!(rs_lcc, A_LCC^2 * γ, c2_lcc, B_LCC^2 * ω, cca2_lcc)
    add_rate!(rs_lcc, A_LCC^3 * γ, c3_lcc, B_LCC^3 * ω, cca3_lcc)
    add_rate!(rs_lcc, A_LCC^4 * γ, c4_lcc, B_LCC^4 * ω, cca4_lcc)
    lcceqs = [ D(x) ~ rs_lcc[x] for x in (c1_lcc, c2_lcc, c3_lcc, c4_lcc, o_lcc, cca0_lcc, cca1_lcc, cca2_lcc, cca3_lcc, cca4_lcc) ]

    rs_ryr = Dict()
    add_rate!(rs_ryr, KA_P_RYR * ca_ss^4, pc1_ryr, KA_M_RYR, po1_ryr)
    add_rate!(rs_ryr, KB_P_RYR * ca_ss^3, po1_ryr, KB_M_RYR, po2_ryr)
    add_rate!(rs_ryr, KC_P_RYR, po1_ryr, KC_M_RYR, pc2_ryr)
    ryreqs = [D(x) ~ rs_ryr[x] for x in (po1_ryr, po2_ryr, pc2_ryr)]

    eqs = [
        1 ~ c0_lcc + c1_lcc + c2_lcc + c3_lcc + c4_lcc + o_lcc + cca0_lcc + cca1_lcc + cca2_lcc + cca3_lcc + cca4_lcc,
        1 ~ pc1_ryr + po1_ryr + po2_ryr + pc2_ryr,
        α_lcc ~ 0.4kHz * exp((v + 2mV) * inv(10mV)),
        β_lcc ~ 0.05kHz * exp(-(v + 2mV) * inv(13mV)),
        y_inf ~ expit(-(v + 55mV) * inv(7.5mV)) + 0.5 * expit((v - 21mV) * inv(6mV)),
        τ_yca ~ 20ms + 600ms * expit(-(v + 30mV) * inv(9.5mV)),
        D(x_yca) ~ (y_inf - x_yca) / τ_yca,
        ICaMax ~ ghk(P_CA_LCC, vm, 1μM, 0.341 * ca_o, 2) ,
        ICaL ~ 6 * x_yca * o_lcc * ICaMax,
        ICaK ~ hil(I_CA_HALF_LCC, ICaMax) * x_yca * o_lcc * ghk(P_K_LCC, vm, k_i, k_o),
        Jrel ~ R_RYR * (po1_ryr + po2_ryr) * (ca_jsr - ca_ss)
    ]

    eqs_cicr16 = [ryreqs; lcceqs; eqs]

    return (; eqs_cicr16, ICaL, ICaK, Jrel)
end
