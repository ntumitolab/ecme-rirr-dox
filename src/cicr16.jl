# Good old 16-state CICR (Cortassa 2006), including L-type calcium channel (LCC) and Ryanodine receptor (RyR)
"LCC ODE system"
function get_lcc_sys(; ca_ss, ca_o, k_i, k_o, vm, name=:cicrsys)
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
        # Currents
        ICaMax(t)
        ICaK(t)
        ICaL(t)
    end

    v = vm / mV
    ω, f, g = ω_LCC, f_LCC, g_LCC
    α = α_lcc
    β = β_lcc
    α′ = 2α
    β′ = 0.5β
    γ = γ_LCC * ca_ss
    vc0c1 = 4α * c0_lcc - β * c1_lcc
    vc1c2 = 3α * c1_lcc - 2β * c2_lcc
    vc2c3 = 2α * c2_lcc - 3β * c3_lcc
    vc3c4 = α * c3_lcc - 4β * c4_lcc
    vc4o = f * c4_lcc - g * o_lcc
    vca0ca1 = 4α′ * cca0_lcc - β′ * cca1_lcc
    vca1ca2 = 3α′ * cca1_lcc - 2β′ * cca2_lcc
    vca2ca3 = 2α′ * cca2_lcc - 3β′ * cca3_lcc
    vca3ca4 = α′ * cca3_lcc - 4β′ * cca4_lcc
    vc0ca0 = γ * c0_lcc - ω * cca0_lcc
    vc1ca1 = 2γ * c1_lcc - ω / 2 * cca1_lcc
    vc2ca2 = 4γ * c2_lcc - ω / 4 * cca2_lcc
    vc3ca3 = 8γ * c3_lcc - ω / 8 * cca3_lcc
    vc4ca4 = 16γ * c4_lcc - ω / 16 * cca4_lcc

    eqs = [
        α_lcc ~ 0.4kHz * exp((v + 2mV) * inv(10mV)),
        β_lcc ~ 0.05kHz * exp(-(v + 2mV) * inv(13mV)),
        1 ~ c0_lcc + c1_lcc + c2_lcc + c3_lcc + c4_lcc + o_lcc + cca0_lcc + cca1_lcc + cca2_lcc + cca3_lcc + cca4_lcc,
        # D(c0_lcc) ~ -vc0c1 - vc0ca0,
        D(c1_lcc) ~ vc0c1 - vc1c2 - vc1ca1,
        D(c2_lcc) ~ vc1c2 - vc2c3 - vc2ca2,
        D(c3_lcc) ~ vc2c3 - vc3c4 - vc3ca3,
        D(c4_lcc) ~ vc3c4 - vc4o - vc4ca4,
        D(o_lcc) ~ vc4o,
        D(cca0_lcc) ~ vc0ca0 - vca0ca1,
        D(cca1_lcc) ~ vc1ca1 + vca0ca1 - vca1ca2,
        D(cca2_lcc) ~ vc2ca2 + vca1ca2 - vca2ca3,
        D(cca3_lcc) ~ vc3ca3 + vca2ca3 - vca3ca4,
        D(cca4_lcc) ~ vca3ca4 + vc4ca4,
        y_inf ~ expit(-(v + 55mV) * inv(7.5mV)) + 0.5 * expit((v - 21mV) * inv(6mV)),
        τ_yca ~ 20ms + 600ms * expit(-(v + 30mV) * inv(9.5mV)),
        D(x_yca) ~ (y_inf - x_yca) / τ_yca,
        ICaMax ~ ghk(P_CA_LCC, vm, 1μM, 0.341 * ca_o, 2) ,
        ICaL ~ 6 * x_yca * o_lcc * ICaMax,
        ICaK ~ hil(I_CA_HALF_LCC, ICaMax) * x_yca * o_lcc * ghk(P_K_LCC, vm, k_i, k_o)
    ]
    return System(eqs, t; name)
end

"Ryanodine receptor (RyR)"
function get_ryr_sys(; ca_jsr, ca_ss, name=:ryrsys)
    @parameters begin
        R_RYR = 3.6kHz
        KA_P_RYR = 1.125E10 / (mM^4 * ms)
        KA_M_RYR = 0.576kHz
        KB_P_RYR = 4.05E6 / (mM^3 * ms)
        KB_M_RYR = 1.93kHz
        KC_P_RYR = 0.1kHz
        KC_M_RYR = 8E-4kHz
    end

    @variables begin
        po1_ryr(t) = 0
        po2_ryr(t) = 0
        pc1_ryr(t) ## Conserved
        pc2_ryr(t) = 0
        Jrel(t)
    end

    vc1o1 = KA_P_RYR * ca_ss^4 * pc1_ryr - KA_M_RYR * po1_ryr
    vo1o2 = KB_P_RYR * ca_ss^3 * po1_ryr - KB_M_RYR * po2_ryr
    vo1c2 = KC_P_RYR * po1_ryr - KC_M_RYR * pc2_ryr

    eqs = [
        1 ~ pc1_ryr + po1_ryr + po2_ryr + pc2_ryr,
        D(po1_ryr) ~ vc1o1 - vo1o2 - vo1c2,
        D(po2_ryr) ~ vo1o2,
        D(pc2_ryr) ~ vo1c2,
        Jrel ~ R_RYR * (po1_ryr + po2_ryr) * (ca_jsr - ca_ss)
    ]
    return System(eqs, t; name)
end
