# Good old CICR, including L-type calcium channel (LCC) and Ryanodine receptor (RyR)
"LCC ODE system"
function get_lcc_sys(ca_ss, ca_o, k_i, k_o, vm; name=:lccsys)
    @parameters begin
        A_LCC = 2
        B_LCC = 2
        ω_LCC = 0.01/ms
        f_LCC = 0.3/ms
        g_LCC = 2/ms
        γ_LCC = 0.1875/(ms * mM)
        P_CA_LCC = 1.24E-3cm * Hz
        P_K_LCC = 1.11E-11cm * Hz
        I_CA_HALF_LCC = -0.4583μA / cm²
    end

    @variables begin
        # State transition rates
        α_lcc(t)
        β_lcc(t)
        # Closed LCC state
        c0_lcc(t) # = 0.9991
        c1_lcc(t) = 8.175e-5
        c2_lcc(t) = 2.508e-9
        c3_lcc(t) = 3.421e-14
        c4_lcc(t) = 1.749e-20
        # Opened LCC state
        o_lcc(t) = 2.624e-20
        # Ca-inhibited LCC state
        cca0_lcc(t) = 1.1328e-3
        cca1_lcc(t) = 1.7591e-8
        cca2_lcc(t) = 6.3826e-13
        cca3_lcc(t) = 1.815e-15
        cca4_lcc(t) = 0.0
        x_yca(t) = 0.9479 # Voltage-gated LCC
        y_inf(t)
        τ_yca(t)
        # Currents
        ICaMax(t)
        ICaK(t)
        ICaL(t)
    end

    v = vm / mV
    A, B = A_LCC, B_LCC
    ω, f, g = ω_LCC, f_LCC, g_LCC
    α = α_lcc
    β = β_lcc
    α′ = α * A
    β′ = β / B
    γ = γ_LCC * ca_ss
    v01 = 4 * α * c0_lcc - β * c1_lcc
    v12 = 3 * α * c1_lcc - 2 * β * c2_lcc
    v23 = 2 * α * c2_lcc - 3 * β * c3_lcc
    v34 = α * c3_lcc - 4 * β * c4_lcc
    v45 = f * c4_lcc - g * o_lcc
    v67 = 4 * α′ * cca0_lcc - β′ * cca1_lcc
    v78 = 3 * α′ * cca1_lcc - 2 * β′ * cca2_lcc
    v89 = 2 * α′ * cca2_lcc - 3 * β′ * cca3_lcc
    v910 = α′ * cca3_lcc - 4 * β′ * cca4_lcc
    v06 = γ * c0_lcc - ω * cca0_lcc
    v17 = A * γ * c1_lcc - ω / B * cca1_lcc
    v28 = (A^2) * γ * c2_lcc - ω / (B^2) * cca2_lcc
    v39 = (A^3) * γ * c3_lcc - ω / (B^3) * cca3_lcc
    v410 = (A^4) * γ * c4_lcc - ω / (B^4) * cca4_lcc

    eqs = [
        α_lcc ~ 0.4kHz * exp((v + 2) / 10),
        β_lcc ~ 0.05kHz * exp(-(v + 2) / 13),
        1 ~ c0_lcc + c1_lcc + c2_lcc + c3_lcc + c4_lcc + o_lcc + cca0_lcc + cca1_lcc + cca2_lcc + cca3_lcc + cca4_lcc,
        # D(c0_lcc) ~ -v01 - v06,
        D(c1_lcc) ~ v01 - v12 - v17,
        D(c2_lcc) ~ v12 - v23 - v28,
        D(c3_lcc) ~ v23 - v34 - v39,
        D(c4_lcc) ~ v34 - v45 - v410,
        D(o_lcc) ~ v45,
        D(cca0_lcc) ~ v06 - v67,
        D(cca1_lcc) ~ v17 + v67 - v78,
        D(cca2_lcc) ~ v28 + v78 - v89,
        D(cca3_lcc) ~ v39 + v89 - v910,
        D(cca4_lcc) ~ v910 + v410,
        y_inf ~ expit(-(v + 55) / 7.5) + 0.5 * expit((v - 21) / 6),
        τ_yca ~ 20ms + 600ms * expit(-(v + 30) / 9.5),
        D(x_yca) ~ (y_inf - x_yca) / τ_yca,
        ICaMax ~ ghk(P_CA_LCC, vm, 1μM, 0.341 * ca_o, 2),
        ICaL ~ 6 * x_yca * o_lcc * ICaMax,
        ICaK ~ hil(I_CA_HALF_LCC, ICaMax) * x_yca * o_lcc * ghk(P_K_LCC, vm, k_i, k_o)
    ]
    return ODESystem(eqs, t; name)
end

"Ryanodine receptor (RyR)"
function get_ryr_sys(ca_jsr, ca_ss; name=:ryrsys)
    @parameters begin
        R_RYR = 3.6 / ms
        KA_P_RYR = 1.125E10 / (mM^4 * ms)
        KA_M_RYR = 0.576 / ms
        KB_P_RYR = 4.05E6 / (mM^3 * ms)
        KB_M_RYR = 1.93 / ms
        KC_P_RYR = 0.1 / ms
        KC_M_RYR = 8E-4 / ms
    end

    @variables begin
        po1_ryr(t) = 8.309e-5
        po2_ryr(t) = 7.12e-11
        pc1_ryr(t) # = 0.7528
        pc2_ryr(t) = 0.2471
        Jrel(t)
    end

    vo1c1 = KA_M_RYR * po1_ryr - KA_P_RYR * NaNMath.pow(ca_ss, 4) * pc1_ryr
    vo1o2 = KB_P_RYR * NaNMath.pow(ca_ss, 3) * po1_ryr - KB_M_RYR * po2_ryr
    vo1c2 = KC_P_RYR * po1_ryr - KC_M_RYR * pc2_ryr

    eqs = [
        1 ~ pc1_ryr + po1_ryr + po2_ryr + pc2_ryr,
        D(po1_ryr) ~ -vo1c1 - vo1o2 - vo1c2,
        D(po2_ryr) ~ vo1o2,
        D(pc2_ryr) ~ vo1c2,
        Jrel ~ R_RYR * (po1_ryr + po2_ryr) * (ca_jsr - ca_ss)
    ]
    return ODESystem(eqs, t; name)
end
