"""
16-state CICR from Cortassa et al. (2006), including L-type calcium channels (LCC) and Ryanodine receptors (RyR)
"""
get_default_cicr16_params() = ComponentArray(;
    R_XFER= inv(0.5747ms),       # Diffusion rate between subspace and cytosol
    ## L-type calcium channels
    A_LCC = 2,
    B_LCC = 1/2,
    ω_LCC = 0.01kHz,
    f_LCC = 0.3kHz,
    g_LCC = 2kHz,
    γ_LCC = 0.1875 / (ms * mM),
    P_CA_LCC = 1E-3cm * Hz, # 1.24E-3cm * Hz
    P_K_LCC = 1.11E-11cm * Hz,
    I_CA_HALF_LCC = -0.4583μAcm⁻²,
    ## Ryanodine receptors
    R_RYR = 3.6kHz,
    KA_P_RYR = 1.125E10 / (mM^4 * ms),
    KA_M_RYR = 0.576kHz,
    KB_P_RYR = 4.05E6 / (mM^3 * ms),
    KB_M_RYR = 1.93kHz,
    KC_P_RYR = 0.1kHz,
    KC_M_RYR = 8E-4kHz,
)

_c0_lcc(c1_lcc, c2_lcc, c3_lcc, c4_lcc, o_lcc, cca0_lcc, cca1_lcc, cca2_lcc, cca3_lcc, cca4_lcc) = 1 - (c1_lcc + c2_lcc + c3_lcc + c4_lcc + o_lcc + cca0_lcc + cca1_lcc + cca2_lcc + cca3_lcc + cca4_lcc)
_c0_lcc(u) = _c0_lcc(u.c1_lcc, u.c2_lcc, u.c3_lcc, u.c4_lcc, u.o_lcc, u.cca0_lcc, u.cca1_lcc, u.cca2_lcc, u.cca3_lcc, u.cca4_lcc)

function cicr16_rates!(D, u, p, t)
    @unpack A_LCC, B_LCC, γ_LCC, ω_LCC, f_LCC, g_LCC = p
    @unpack vm, ca_ss, c1_lcc, c2_lcc, c3_lcc, c4_lcc, o_lcc, cca0_lcc, cca1_lcc, cca2_lcc, cca3_lcc, cca4_lcc = u

    ## Calcium activation / inactivation rates
    ω, f, g = ω_LCC, f_LCC, g_LCC
    α = 0.4kHz * exp((vm + 2mV) * inv(10mV))
    β = 0.05kHz * exp(-(vm + 2mV) * inv(13mV))
    α′ = A_LCC * α
    β′ = B_LCC * β
    γ = γ_LCC * ca_ss
    c0_lcc = _c0_lcc(u)
    vc0c1 = 4α * c0_lcc - β * c1_lcc
    vc1c2 = 3α * c1_lcc - 2β * c2_lcc
    vc2c3 = 2α * c2_lcc - 3β * c3_lcc
    vc3c4 = α * c3_lcc - 4β * c4_lcc
    vc4o = f * c4_lcc - g * o_lcc
    vcca0cca1 = 4α′ * cca0_lcc - β′ * cca1_lcc
    vcca1cca2 = 3α′ * cca1_lcc - 2β′ * cca2_lcc
    vcca2cca3 = 2α′ * cca2_lcc - 3β′ * cca3_lcc
    vcca3cca4 = α′ * cca3_lcc - 4β′ * cca4_lcc
    vc0cca0 = γ * c0_lcc - ω * cca0_lcc
    vc1cca1 = A_LCC * γ * c1_lcc - B_LCC * ω * cca1_lcc
    vc2cca2 = A_LCC^2 * γ * c2_lcc - B_LCC^2 * ω * cca2_lcc
    vc3cca3 = A_LCC^3 * γ * c3_lcc - B_LCC^3 * ω * cca3_lcc
    vc4cca4 = A_LCC^4 * γ * c4_lcc - B_LCC^4 * ω * cca4_lcc
    D.c1_lcc = vc0c1 - vc1c2 - vc1cca1
    D.c2_lcc = vc1c2 - vc2c3 - vc2cca2
    D.c3_lcc = vc2c3 - vc3c4 - vc3cca3
    D.c4_lcc = vc3c4 - vc4o - vc4cca4
    D.o_lcc = vc4o
    D.cca0_lcc = vc0cca0 - vcca0cca1
    D.cca1_lcc = vcca0cca1 - vcca1cca2 + vc1cca1
    D.cca2_lcc = vcca1cca2 - vcca2cca3 + vc2cca2
    D.cca3_lcc = vcca2cca3 - vcca3cca4 + vc3cca3
    D.cca4_lcc = vcca3cca4 + vc4cca4


end

function get_cicr16_eqs(; vm, ca_ss, ca_i, ca_jsr, ca_nsr, k_i, ca_o, k_o)
    @parameters begin
        R_TR = inv(9.09ms)          # Diffusion rate between JSR and NSR
        R_XFER= inv(0.5747ms)       # Diffusion rate between subspace and cytosol
        ## L-type calcium channels
        A_LCC = 2
        B_LCC = 1/2
        ω_LCC = 0.01kHz
        f_LCC = 0.3kHz
        g_LCC = 2kHz
        γ_LCC = 0.1875 / (ms * mM)
        P_CA_LCC = 1E-3cm * Hz # 1.24E-3cm * Hz
        P_K_LCC = 1.11E-11cm * Hz
        I_CA_HALF_LCC = -0.4583μAcm⁻²
        ## Ryanodine receptors
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
        c0_lcc(t) # Conserved
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
        Jtr(t)
        Jxfer(t)
    end

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
        α_lcc ~ 0.4kHz * exp((vm + 2mV) * inv(10mV)),
        β_lcc ~ 0.05kHz * exp(-(vm + 2mV) * inv(13mV)),
        y_inf ~ expit(-(vm + 55mV) * inv(7.5mV)) + 0.5 * expit((vm - 21mV) * inv(6mV)),
        τ_yca ~ 20ms + 600ms * expit(-(vm + 30mV) * inv(9.5mV)),
        D(x_yca) ~ (y_inf - x_yca) / τ_yca,
        ICaMax ~ ghk(P_CA_LCC, vm, 1μM, 0.341 * ca_o, 2) ,
        ICaL ~ 6 * x_yca * o_lcc * ICaMax,
        ICaK ~ hil(I_CA_HALF_LCC, ICaMax) * x_yca * o_lcc * ghk(P_K_LCC, vm, k_i, k_o),
        Jrel ~ R_RYR * (po1_ryr + po2_ryr) * (ca_jsr - ca_ss),
        Jtr ~ R_TR * (ca_nsr - ca_jsr),
        Jxfer ~ R_XFER * (ca_ss - ca_i),
    ]

    eqs_cicr16 = [ryreqs; lcceqs; eqs]

    return (; eqs_cicr16, ICaL, ICaK, Jrel, Jtr, Jxfer)
end

function get_cicr16_sys(; vm, ca_ss, ca_i, ca_jsr, ca_nsr, k_i, ca_o, k_o, name=:cicr16sys)
    @unpack eqs_cicr16 = get_cicr16_eqs(; vm, ca_ss, ca_i, ca_jsr, ca_nsr, k_i, ca_o, k_o)
    return System(eqs_cicr16, t; name)
end
